"predict.obliqueRF" <-
function(object, newdata, type="response", proximity=F, ...)
{
    if (!inherits(object, "obliqueRF"))
        stop("object not of class obliqueRF")
    if (is.null(object$trees)) stop("No trees in the object")
    if ( object$type != "classification" )
        stop("right now predict supports only (binary) classification")

	typeOK=switch(type,response="ok",prob="ok",votes="ok","notok");
	if(typeOK=="notok")
		stop("specified type is not valid, choose one of \"prob\", \"votes\" or \"response\"; aborting...");
		
    ntree=object$ntree;
    num_classes=object$num_classes
    num_samples=nrow(newdata);
    newdata_class_votes<-matrix(0, num_classes, num_samples);

    # data frames are just no fun
    if(!is.matrix(newdata))
      newdata<-as.matrix(newdata);
    
    if(proximity)
    	proxM=matrix(NA,num_samples,ntree);
    
    #loop over all trees
    for(currTree in 1:ntree)
    {
      #for some very rough profiling
      #cat(currTree,"th tree of",ntree,date(),"\n");
		tree=object$trees[[currTree]];
		# convert the node_class from 0/1 to 1/2
		node_classIdx=tree$node_class+1;

		for(M in 1:num_samples)
		{
			curr_node=1;
			#go through all nodes of current tree
			for(i in 1:(tree$total_nr_nodes))
			{
				if(tree$node_status[curr_node]==-1) #==terminal node
				{
					# increment class count by 1
					newdata_class_votes[node_classIdx[curr_node],M]=
						newdata_class_votes[node_classIdx[curr_node],M]+1;
						
				    if(proximity)
						proxM[M,currTree]<-curr_node

					break;
				}
				dset<-newdata[M, tree$best_vars[,curr_node]];
				score=sum(tree$training_stats[,curr_node]*dset);
				if(score<tree$best_split_value[curr_node])
				{
					curr_node=tree$treemap[1, curr_node]
				}
				else
				{
					curr_node=tree$treemap[2, curr_node];
				}
			}
		}
    }
    # count the total number of votes; it's the same for all
    #  rows - so why not just take the first one...
    #totalVotes=sum(newdata_class_votes[,1]);
    # should be the same as the number of trees...
    totalVotes=ntree;
    
    # normalise the votes for the second class by the total
    #  number of votes; output can be tested for >.5 or
    #  something similar
    pred=c();
    if(type=="response")
		pred=object$class_names[1+((newdata_class_votes[2,]/totalVotes)>0.5)*1];
    if(type=="prob"){
		pred=t(newdata_class_votes/totalVotes);
        colnames(pred)<-object$class_names
    }
    if(type=="votes"){
		pred=t(newdata_class_votes);
        colnames(pred)<-object$class_names
    }
	
	if(!proximity)
	{
		res<-pred;
	}
	else
	{
		node_counts<-matrix(0,num_samples,num_samples)
		#check for every tree the node numbers and check what samples
		#have common ones
		for(op in 1:ntree){
			tabl<-table(proxM[,op])
			#get different final node numbers for current tree
			mults<-as.numeric(names(tabl)[tabl>1])
			#for every node number increase the count of the respective samples
			for(mi in mults){
				wh<-which(proxM[,op]==mi)
				node_counts[wh, wh]<-node_counts[wh, wh] +1
			}
		}
		#normalize by the number of trees
		node_counts<-node_counts/ntree;
		
		res<-list( pred=pred, proximity=node_counts);
	}
		
	return(res);
}
