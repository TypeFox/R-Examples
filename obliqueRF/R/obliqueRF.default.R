

##--------------------------------------
## oblique random forest 
##--------------------------------------

obliqueRF.default<- function(
                    x,
                    y,
                    x.test=NULL,
                    y.test=NULL,
                    mtry=NULL,
                    ntree=100,
                    training_method="ridge",                    
                    bImportance=F,
                    bProximity = F,
                    verbose=F,
                    ...
                    )
{

    num_split_vars <- mtry

    # always check out-of-bag error
    do_oob<-TRUE

    if(verbose)
      cat("obliqueRF: WARNING: status messages enabled; set \"verbose\" to false to disable\n");

	# --- importance ---

	impStruct<-list(bImportance=bImportance,importance_pvalue=0,importance_weight="");

    if(bImportance){
        impStruct$importance_pvalue <- 0.05 # p value
        impStruct$importance_weight <- "none"
        do_oob<-FALSE  
        
        if(bProximity){
            cat("obliqueRF: WARNING: importance and proximity cannot be calculated at the same time. Will set \'bProximity = FALSE\' and continue.\n")
            bProximity = FALSE
        }

        if(!(training_method=="log")){
            cat("obliqueRF: WARNING: importance requires training_method=\"log\". Adjusting training_method...\n");
            training_method="log"    
        }

    }
    
    # --- check training method selection ---

    node_train_method=NULL;    
    if(training_method=="ridge"){        
        node_train_method=ridge_svd_wrapper
        cat("obliqueRF: using fast ridge regression with SVD as node model,\n does not scale data at the node,\n tests regularization parameter \"lambda\" for 10^(c(-5:5)).\n In case of problems with the SVD variant try training_method=\"ridge_slow\" for separate ridge regressions (slower).\n")
    }
    if(training_method=="ridge_slow"){        
        node_train_method=ridge_slow_wrapper
        cat("obliqueRF: using slow ridge regression as node model,\n does not scale data at the node,\n tests regularization parameter \"lambda\" for 10^(c(-5:5)) with separate explicit ridge regressions.\n")
    }
    if(training_method=="log"){
        node_train_method=log_wrapper
        cat("obliqueRF: using logistic regression as node model.\n")        
    }
    if(training_method=="svm"){
        node_train_method=svm_wrapper        
        cat("obliqueRF: using linear SVM as node model,\n scales data at the node,\n tests regularization parameter \"cost\" for 10^(c(-3,-1,1)).\n")        
    }
    if(training_method=="pls"){
        node_train_method=pls_wrapper        
        cat("obliqueRF: using partial least squares regression as node model,\n does not scale data at the node,\n optimizes regularization parameter \"ncomp\" (number of components) for 1:min(20, mtry).\n")
    }        
    if(training_method=="rnd"){        
        node_train_method=rnd_wrapper
        cat("obliqueRF: using a random hyperplane as node model,\n scales data to unit variance at the node,\n coefficients of random subspace projection are chosen from a normal distribution with unit variance.\n")
    }
	if( is.null(node_train_method) )
	{
		stop("obliqueRF: the training method is not valid, choose one of the following:\n  ridge, ridge_slow, log, svm, pls or rnd.\n  See documentation for more details. Aborting");
	}

    node_pred_method=predict

    # generic impurity support
    node_impurity_method=gini_impurity_wrapper;


    # ------------- DATA -------------

    # make sure we have valid data
    if(
      any(is.nan(as.vector(x))) |
      any(is.infinite(as.vector(x))) |
      any(is.nan(as.vector(y))) |
      any(is.infinite(as.vector(y)))
      )
    stop("obliqueRF: Training data contains NaNs or Inf values. Please correct and try again. Aborting...");

    if(is.null(y.test)){
        cat("obliqueRF: no test set defined. Will use training data.\n")
        x.test<-x
        y.test<-y    
    }
    if(
      sum(is.nan(as.vector(x.test)))>0 |
      sum(is.infinite(as.vector(x.test)))>0 |
      sum(is.nan(as.vector(y.test)))>0 |
      sum(is.infinite(as.vector(y.test)))>0
      )
      stop("obliqueRF: Test data contains NaNs or Inf values. Please correct and try again. Aborting...");

    # labels: transform factors etc. to numerical values in [1,0]
    clnames <- sort(unique(paste(y)))
    if(!all(y==0|y==1)){
        
         yout<-rep(0,length(y))
         ytestout<-rep(0,length(y.test))
         
         yout[y==clnames[1]]<-0
         ytestout[y.test==clnames[1]]<-0     
         yout[y==clnames[2]]<-1
         ytestout[y.test==clnames[2]]<-1     
        
         y<-yout
         y.test<-ytestout
    }    

    # get a data.frame from x, y
    num_vars=ncol(x);
    num_samples=nrow(x);    

    if(length(y)!=num_samples)
        stop("obliqueRF: WARNING: number of samples and length of y does not match.")

    test_set<-cbind(as.matrix(x.test),y.test)
    
    d<-cbind(as.matrix(x),y)
    cl_id=num_vars+1;
    cl<- y

    # number of variables to split on at each node
    if(is.null(num_split_vars))
        num_split_vars=max(as.integer(sqrt(num_vars)),2); # default
    if(num_split_vars> num_vars){
        stop('obliqueRF: WARNING: mtry > number of variable.')
    }

    # check number of variables
    num_classes<-length(unique(cl));
    if((num_classes<2) || min(cl))
        stop("obliqueRF: WARNING: wrong number of classes. Aborting...")

    # ------------- VARIABLES -----------------

    # allocate variables
      nr_nodes=2*num_samples+1; #total number of possible nodes
      nd_min_size=1; #minimum size of a node
      oob_glob=rep(0, num_samples); #indicates how often a sample was oob

      oob_class_votes <- matrix(0, num_classes, nrow(d));
      test_set_class_votes<-matrix(0, num_classes, nrow(test_set));

      if(bImportance){
		  imp<-rep(0,ncol(x))
		  imp_class_votes<-matrix(0, num_classes, nrow(d));
      }

      oob_prox <- matrix(NA, nrow(d), ntree)

    # ------------- TRAIN -------------

     err<-c()
     err_0<-c()
     err_1<-c()    
    
    # main loop
    trees=vector("list",ntree);
	if(verbose)
        cat("obliqueRF: starting.",date(),"\n")
    
    for(curr_num_tree in 1:ntree)
        {   
			if(verbose)
			  if((curr_num_tree%%10)==1) cat("", curr_num_tree,"out of",ntree,"trees so far.",date(),"\n");

            class_n=rep(0, num_classes); #vector containing the class distribution in sample

            # sample with replacement
            k<-sample(1:num_samples, num_samples, replace=TRUE);

            # "bag" statistic, how often is a sample in the bag?
            sample_pop=as.numeric(table(factor(k,levels=1:num_samples)));
            # count the class distribution in sample (weights here?)
            class_n=as.numeric(table(cl[k]));

            dframe<-d
            b_sample=sample_pop; # sampling frequency and index vector for oob
            nr_diff_cases=length(sample_pop>0); #nr of different (unique) cases in the bootstrap sample
            n_case=1:num_samples;   #translation to original case numbers (needed for class population, etc.)

            # ---------------------- BUILDTREE ------------------------
            # build the tree, i.e. nr of columns in dframe (-1 for cl)

            tree<-buildtree(
                dframe,
                b_sample,
                class_n,
                n_case,
                nr_nodes,
                nd_min_size,
                nr_diff_cases,
                num_vars,
                num_split_vars,
                num_classes,
                node_train_method,
                node_pred_method,
				node_impurity_method,
                num_samples,
                impStruct
                );
            trees[[curr_num_tree]]=tree;

            # ---------------------- PREDICT ------------------------

			# convert the node_class from 0/1 to 1/2
			node_classIdx=tree$node_class+1;

            # -----------------
            # oob error
            if(do_oob){
				#global: count the oob's
				oob_glob=oob_glob+(sample_pop==0)*1;

				for(M in which(sample_pop==0)) # only oob samples
				{
					curr_node=1;
					for(i in 1:(tree$total_nr_nodes))
					{
						if(tree$node_status[curr_node]==-1) #==terminal node
						{
							# increment class count by 1
							oob_class_votes[node_classIdx[curr_node],M]=
							  oob_class_votes[node_classIdx[curr_node],M]+1;
							  
							oob_prox[M,curr_num_tree]<-curr_node
							break;
						}
						dset<-d[M, tree$best_vars[,curr_node]]; # *NOT* test_set
						score=sum(tree$training_stats[,curr_node]*dset);
						if(score<tree$best_split_value[curr_node])
							curr_node=tree$treemap[1, curr_node]
						else
							curr_node=tree$treemap[2, curr_node];
					}
				}
            }
            # -----------------
            # importance

            if(bImportance){
                for(pi in 1:ncol(tree$signif))
                 imp[tree$best_vars[,pi]]<-imp[tree$best_vars[,pi]] + (tree$significance[,pi]) 
                 #plot(imp,t="l", xlab="Features", ylab="Importance", main=curr_num_tree)
            }

            # -----------------
            # test set?
            
                for(M in 1:nrow(test_set))
                {
                    curr_node=1;
                    for(i in 1:(tree$total_nr_nodes))
                    {
                        if(tree$node_status[curr_node]==-1) #==terminal node
                        {
						    # increment class count by 1
                            test_set_class_votes[node_classIdx[curr_node],M]=
						      test_set_class_votes[node_classIdx[curr_node],M]+1;
                            break;
                        }
                        dset<-test_set[M, tree$best_vars[,curr_node]];
                        score=sum(tree$training_stats[,curr_node]*dset);
                        if(score<tree$best_split_value[curr_node])
                            curr_node=tree$treemap[1, curr_node]
                        else
                            curr_node=tree$treemap[2, curr_node];
                    }
                }

            # -------------------- EVAL --------------------

            v<- 1*(test_set_class_votes[2,]/sum(test_set_class_votes[,1])>.5)
            m<-table(factor(y.test,levels=c(0,1)),factor(v,levels=c(0,1)))
            err_0[curr_num_tree]<- m[1,2]/sum(m[1,])
            err_1[curr_num_tree]<- m[2,1]/sum(m[2,])
            err[curr_num_tree]<- (m[1,2]+m[2,1])/sum(m)


        
     } # loop over trees

    test.votes<-test_set_class_votes[2,]/sum(test_set_class_votes[,1])

    # calc proximity
    if(bProximity){
		oob_counts<-matrix(0,nrow(d),nrow(d)) # diag(oob_counts) is oob_glob
		#count how often each sample was oob
		for(op in 1:ncol(oob_prox)){                
			wh<-which(!is.na(oob_prox[,op]))
			oob_counts[wh, wh]<-oob_counts[wh, wh] +1
		}

		node_counts<-matrix(0,nrow(d),nrow(d))
		for(op in 1:ncol(oob_prox)){
			tabl<-table(oob_prox[,op])
			mults<-as.numeric(names(tabl)[tabl>1])
			for(mi in mults){
			wh<-which(oob_prox[,op]==mi)
			node_counts[wh, wh]<-node_counts[wh, wh] +1
			}
		} 
		dst<-1 - node_counts/oob_counts
		diag(dst)<-0
    }



    oob.votes<-NULL; thr<-NULL
    if(do_oob){
		oob.votes<-oob_class_votes[2,]/(oob_class_votes[1,]+oob_class_votes[2,])
		pred <- prediction(oob.votes, y)
		perf <- performance(pred,"acc")
		thr<-perf@"x.values"[[1]][which.max(perf@"y.values"[[1]])]
    }

    cl <- match.call()
    cl[[1]] <- as.name("obliqueRF")

    res<-list(
    	call=cl,
    	type = "classification",
        class_names = clnames,
    	errs=list(err_0 = err_0, err_1 = err_1, err = err),
    	pred=list(  pred.test = 1*(test.votes>.5),
    		    votes.test = test.votes,
    		    y.test = y.test,
    		    votes.oob = oob.votes,
    		    y = y),
    	lab=paste('oRF - node:', tree$lab),
    	ntree = ntree,
    	mtry = num_split_vars,
    	importance = if(bImportance) imp else NULL,
    	proximity = if(bProximity) dst else NULL,
    	num_classes = num_classes,
    	trees = trees
    )

    class(res) <- "obliqueRF"
    return(res);
}


##--------------------------------------------
## buildtree
##--------------------------------------------


buildtree <- function(
    dframe,
    b_sample,
    class_n,
    n_case,
    nr_nodes,
    nd_min_size,
    nr_diff_cases,
    num_vars,
    num_split_vars,
    num_classes,
    node_train_method,
    node_pred_method,
    node_impurity_method,
    num_samples,
    impStruct
    )
{
    class_distr_nodes=matrix(0, num_classes, nr_nodes); #class distribution per node
    node_status=rep(0, nr_nodes); #indicates what type a node is: 1: parent node, 2: has not been processed yet, -1: terminal node
    node_start=rep(0, nr_nodes);  #start index for a node in the data matrix
    node_length=rep(0, nr_nodes); #nr of cases belonging to that node, starting from node_start

    var_imp=rep(0, num_vars);

    best_vars=matrix(0, num_split_vars, nr_nodes); #holds the split variables
    treemap<-matrix(0,2,nr_nodes); #holds the indexes of the child nodes
    dgini=rep(0, nr_nodes);
    dgini_perc=rep(0, nr_nodes);
    significance=matrix(0, num_split_vars, nr_nodes);
    best_split_values=rep(0,nr_nodes);
    ridge_values=rep(NA,nr_nodes);
    training_stats=matrix(0, num_split_vars,nr_nodes);

    ncur=1; #index of the node we can write to
    node_status[1]=2; #first node has not been split yet...
    node_start[1]=1;
    node_length[1]=nr_diff_cases;

    #copy the class-numbers
    class_distr_nodes[,1]=class_n;

    n_temp_case=n_case;
    #iterate over the node_status array and check if there are unprocessed nodes;
    #split if necessary until terminal nodes are reached
    nr_parent_nds=0;
    curr_node=0;
    lab<-" "
    
    for(curr_node in 1:nr_nodes){

        if(node_status[curr_node]==0) break;
        if(node_status[curr_node]!=2) next; #the current node is either a parent or a terminal node
        
        #start and end index of the cases in the data matrix
        nd_start=node_start[curr_node];
        nd_length=node_length[curr_node];
        nd_end=nd_start+nd_length-1;

        #get the class distribution for the current node
        class_n=class_distr_nodes[,curr_node];

        res<-try(findbestsplit(
           dframe=dframe[n_case[nd_start:nd_end],],
            b_sample=b_sample[n_case[nd_start:nd_end]],
            class_n,
            nd_start,
            nd_length,
            num_vars,
            num_split_vars,
            num_classes,
            node_train_method,
            node_pred_method,
		    node_impurity_method,
            impStruct
            ),TRUE);

        if((res$jstat==1)|(class(res)=="try-error"))
        {
            node_status[curr_node]=-1;
            next;
        }

        lab<-res$lab

        best_vars[,curr_node]=res$variables;
        dgini[curr_node]=res$decsplit;
        dgini_perc[curr_node]=res$decsplit_perc;
        significance[,curr_node]=res$signif;
        best_split_values[curr_node]=res$best_split_value;
        ridge_values[curr_node]=res$ridge_value;
        training_stats[,curr_node]=res$training_stats;
        var_imp[res$variables]=var_imp[res$variables]+abs(res$training_stats)*((res$decsplit));#(nd_length/nr_diff_cases))

        #count the parent nodes
        nr_parent_nds=nr_parent_nds+1;

        #save size, etc. for the two child nodes
        node_length[ncur+1]=res$isplit;
        node_length[ncur+2]=node_length[curr_node]-res$isplit;
        node_start[ncur+1]=nd_start;
        node_start[ncur+2]=nd_start+res$isplit;

        #update n_case
        true_left_ind=(nd_start-1)+((1:node_length[curr_node])[res$split_vect==1]);
        true_right_ind=(nd_start-1)+((1:node_length[curr_node])[res$split_vect==0]);
        n_temp_case[nd_start:(nd_start+res$isplit-1)]=n_case[true_left_ind];
        n_temp_case[(nd_start+res$isplit):nd_end]=n_case[true_right_ind];
        n_case[nd_start:nd_end]=n_temp_case[nd_start:nd_end];

        #class population/distribution  # sum over b_sample
        for(n in nd_start:(nd_start+res$isplit-1))
        {
            nc=n_case[n];
            j=dframe[nc, num_vars+1]+1;
            class_distr_nodes[j,ncur+1]=class_distr_nodes[j,ncur+1]+b_sample[nc];
        }
        for(n in (nd_start+res$isplit):nd_end)
        {
            nc=n_case[n];
            j=dframe[nc, num_vars+1]+1;
            class_distr_nodes[j,ncur+2]=class_distr_nodes[j,ncur+2]+b_sample[nc];
        }

        #adjust the node status
        node_status[ncur+1]=2;
        node_status[ncur+2]=2;
        if (node_length[ncur+1]<=nd_min_size) node_status[ncur+1]=-1;
        if (node_length[ncur+2]<=nd_min_size) node_status[ncur+2]=-1;

        #if pure node then terminal node
        if(sum(class_distr_nodes[,ncur+1]>0)==1) node_status[ncur+1]=-1;
        if(sum(class_distr_nodes[,ncur+2]>0)==1) node_status[ncur+2]=-1;
        nr_samples=sum(class_distr_nodes[,ncur+1])+sum(class_distr_nodes[,ncur+2]);

        #remember the child nodes
        treemap[1,curr_node]=ncur+1;
        treemap[2,curr_node]=ncur+2;

        node_status[curr_node]=1; #current node has been processed
        ncur=ncur+2;
    } #for(curr_node in 1:nr_nodes)

    #normalise the variable importance
    var_imp=var_imp/nr_parent_nds;
    if(nr_parent_nds>0) var_imp=var_imp/nr_parent_nds;

    total_nr_nodes=curr_node;
    node_status[node_status==2]=-1;

    #get the class for each node
    node_class=apply(class_distr_nodes[,1:total_nr_nodes],2,which.max)-1;

    res2<-list(
        best_vars=best_vars[,1:total_nr_nodes], 
        best_split_values=best_split_values[1:total_nr_nodes], 
        ridge_values=ridge_values[1:total_nr_nodes], 
        training_stats=training_stats[,1:total_nr_nodes], 
        significance=significance[,1:total_nr_nodes],
        node_class=node_class, 
        dgini=dgini[1:total_nr_nodes], 
        dgini_perc=dgini_perc[1:total_nr_nodes], 
        total_nr_nodes=total_nr_nodes, 
        treemap=treemap[,1:total_nr_nodes], 
        node_status=node_status[1:total_nr_nodes], 
        var_imp=var_imp,
        lab=lab
        );

    return(res2);
}


##--------------------------------------------
## findbestsplit
##--------------------------------------------


findbestsplit <- function(
    dframe,
    b_sample,
    class_n,
    nd_start,
    nd_length,
    num_vars,
    num_split_vars,
    num_classes,
    node_train_method,
    node_pred_method,
    node_impurity_method,
    impStruct
    )
{
    #--------

    #draw some variable indices without replacement
    vars<- sample(1:num_vars, num_split_vars);

    #construct the subset (don't forget the labels!)
    temp_matrix<-dframe[,vars];
    cls<-dframe[,num_vars+1];



    #do the training 
    #we can get warnings from the training function - ignore for now
      ow<-options("warn");
	  options(warn=-1);
      z<- try(node_train_method(x=temp_matrix, y=cls, ind = b_sample, impStruct = impStruct),TRUE);
#      z<- try(node_train_method(x=temp_matrix, y=cls, ind = b_sample, impStruct = impStruct),F);

      options(ow);

      lab<-" "
      
      if((class(z)=="try-error") || (sum(is.na(z$coef))>0) || (table(cls)==1) || is.null(z$coef) )
      {
	  jstat=1;
	  res<-list(jstat=jstat);
	  return(res);
      }

    #call some impurity function
    impurity_meas<-node_impurity_method(temp_matrix, cls, b_sample, class_n, num_classes, nd_length, z);

    #if something went wrong we have no valid impurity data
    # the jstat will be tested by the calling function as well
    if( impurity_meas$jstat==1 )
    return(impurity_meas);


    #assign the results
    res<-list(
        jstat=impurity_meas$jstat,
        variables=vars,
        isplit=impurity_meas$isplit,
        split_vect=impurity_meas$split_vect,
        decsplit=impurity_meas$decsplit,
        decsplit_perc=impurity_meas$decsplit_perc,
        best_split_value=impurity_meas$best_split_value,
        ridge_value=z$ridge_value,
        training_stats=z$coef,
        signif=z$signif,
        lab=z$lab
        );
    return(res);
}



##--------------------------------------------
## Impurity wrapper
##--------------------------------------------
gini_impurity_wrapper<-function(dmatrix, cls, b_sample, class_n, num_classes, nd_length, z,...)
{
    #coefficients for gini index
    pno<-0;#numerator
    pdo<-0;#denominator

    pno<-sum(class_n^2);
    pdo<-sum(class_n);

    crit0<-pno/pdo;
    critmax<-(-1.0e20);
    jstat<-0;
    isplit<-0;


    #get the scores
#     pred<-as.matrix(dmatrix)%*%z$coef;
    pred<-(dmatrix)%*%z$coef;
    #sort them
    pred_sort_matching<-sort.list(pred);
    p<-pred[pred_sort_matching];

    uv<-b_sample[pred_sort_matching];
    clv<-cls[pred_sort_matching]+1;

    rldv<-cumsum(uv);
    rrdv<-pdo-rldv;

    cum_class<-array(0,c(num_classes,nd_length,2));
    temp_uv<-numeric(nd_length);
    for(i in 1:num_classes)
    {
        idv<-clv==i;
        temp_uv[idv]<-uv[idv];
        temp_cumsum<-cumsum(temp_uv);
        cum_class[i,,1]<-(temp_cumsum^2);
        cum_class[i,,2]<-((class_n[i]-temp_cumsum)^2);
        temp_uv<-numeric(nd_length);
    }
    rlnv<-apply(cum_class[,,1], 2, sum);
    rrnv<-apply(cum_class[,,2], 2, sum);

    crit<-(rlnv/rldv)+(rrnv/rrdv);
    isplit<-which.max(crit)
    critmax<-crit[isplit];

    if (critmax<(-1.0e10))
    {
        jstat<-1;
        res<-list(jstat=jstat);
        return(res);
    }

    split_vect<-numeric(nd_length);
    #assign 1 to the first group, the 2nd will stay 0
    split_vect[pred_sort_matching[1:isplit]]<-1;
    decsplit<-((critmax-crit0));

    #calculate the offset
    best_split_value<-(p[isplit]+p[isplit+1])/2;
    #assign the results
    res<-list(
	jstat=jstat,
        isplit=isplit,
        split_vect=split_vect,
        decsplit=decsplit,
        best_split_value=best_split_value,
        decsplit_perc=(decsplit/crit0)
        );
    return(res);
}





## -------------------------------------------
## NODE MODELS
## -------------------------------------------


##--------------------------------------------
## Logistic wrapper
## if "IMPORTANCE" is chosen ANOVA is calculated at every split
##--------------------------------------------

# require(stats)

log_wrapper<-function(x=NULL, y=NULL, ind, impStruct,...) {

   ## what samples: get bootstrapped data set 
   wh<-c();
   for(i in 1:max(ind)) wh<-c(wh,which(ind>=i))

   ## what features: remove constant variables
   wp<-which(unlist(lapply(1:ncol(x),function(A) length(unique(x[wh,A]))>1)))

    X<-scale(x[wh,wp])
    z<-try(glm(y[wh]~.,data=as.data.frame(X),family=binomial(link = "logit")),TRUE)

    if(class(z)[1]!="try-error")
    {
        COEF<-coef(z)[-1]/attr(X,"scaled:scale")  

        v_norm=sqrt(COEF%*%COEF);
        z$coef<-rep(0,ncol(x))
        z$coef[wp]<-COEF/v_norm;

	    z$signif<-rep(0,ncol(x))
		if(impStruct$bImportance)
		{
		    if(impStruct$importance_weight == "none")  
		     z$signif[wp]<- c(-1,1)[1+(anova(z, test="Chisq")[-1,5]<impStruct$importance_pvalue)*1]
		    if(impStruct$importance_weight == "length")  
		     z$signif[wp]<- c(-1,1)[1+(anova(z, test="Chisq")[-1,5]<impStruct$importance_pvalue)*1]*length(wh)
		}

        z$lab<-"logit"
        z$ridge_value = 0
    }
    return(z);
}



##--------------------------------------------
## PLS  wrapper 
##--------------------------------------------



pls_wrapper<-function(x=NULL, y=NULL, ind,...) {

       ## what samples: get bootstrapped data set 
       wh<-c();
       for(i in 1:max(ind)) wh<-c(wh,which(ind>=i))
    
       ## what features: remove constant variables
       wp<-which(unlist(lapply(1:ncol(x),function(A) length(unique(x[wh,A]))>1)))

    m_steps<-min(20,length(wp));

    lbd<-c()
    lab<-c(-1,1)[y+1]
    
    X<-scale(x[wh,wp], scale=F)

    obj<-kernelpls.fit(X=as.matrix(X), Y=lab[wh], ncomp =m_steps, stripped = T)

    # find optimal lambda from oob
    Xt<-scale(x[ind==0,wp], center=attr(X,"scaled:center"), scale=F)
    for(l in 1:m_steps){
      if(sum(ind==0)>1)
      lbd[l]<-sum((as.matrix(Xt)%*%coef(obj)[,,l]-as.numeric(obj$Xm%*%coef(obj)[,,l])>0)*1==y[ind==0])
    }
    
    wm<-max(which.max(as.numeric(lbd))[1],2)
    z<-obj

    COEF<-coef(z)[,,wm]

    v_norm=sqrt(COEF%*%COEF);
    z$coef<-rep(0,ncol(x))
    z$coef[wp]<-COEF/v_norm;

    z$signif<-rep(0,ncol(x))
    z$signif[wp]<-0;

    z$ridge_value = wm
    z$lab<-"PLS"
        
    return(z);
}


##--------------------------------------------
## ridge wrapper
##--------------------------------------------


ridge_slow_wrapper<-function(x=NULL, y=NULL, ind,...) {
   ## what samples: get bootstrapped data set 
    # ind contains the sample population;
    # build the bootstrap vector with repetitions
    wh<-c();
    for(i in 1:max(ind)) wh<-c(wh,which(ind>=i))

   ## what features: remove constant variables
    wp<-which(unlist(lapply(1:ncol(x),function(A) length(unique(x[wh,A]))>1)))

    numLambdas=11;
    ridgeValues=c(-5:5);
    lambdas=10^(ridgeValues); #vector of possible lambdas
    
    bestRidgeValue=ridgeValues[1];

    # the current sample and feature subset for the node
    X<-scale(x[wh,wp], scale=F)

    lbd<-rep(0, numLambdas); #vector for number of matches per lambda
    lab<-c(-1,1)[y+1] # convert two class case into fit to -1 and 1

   ## if there are no oobs left there is no need for testing different lambdas
    if(sum(ind==0)==0)
    {
      z<-gen.ridge(x=t(t(X)), y=lab[wh], lambda=lambdas[1]);
    }
    else #test different lambdas
    {
      #save objects
      objs<-vector("list",numLambdas);
      
      # the oobs with the current features
      Xt<-scale(x[ind==0,wp],  center=attr(X,"scaled:center"), scale=F)

      # loop over lambda
      for(l in 1:numLambdas){
	  
	  objs[[l]]<-gen.ridge(x=t(t(X)), y=lab[wh], lambda=lambdas[l])
	  p<-predict(objs[[l]], t(t(Xt)));

	  # find optimal lambda from (local) oob: check how many classes match
	  lbd[l]<-sum((p>0)*1==y[ind==0])
	  
      }
      # find (first) best lambda (with best matching)
      wm<-which.max(as.numeric(lbd))[1]
      # get again the stats for the best lambda
      #z<-gen.ridge(x=t(t(X)), y=lab[wh], lambda=lambdas[wm])
      z<-objs[[wm]];
      bestRidgeValue=ridgeValues[wm];
    }

    COEF<-z$coefficients[,1]

    v_norm=sqrt(COEF%*%COEF);
    z$coef<-rep(0,ncol(x))
    z$coef[wp]<-COEF/v_norm;

    z$signif<-rep(0,ncol(x))
    z$signif[wp]<-0;

    z$ridge_value = bestRidgeValue;
    z$lab<-"ridge"

    return(z);
}


##--------------------------------------------
## ridge wrapper using SVD
##--------------------------------------------


ridge_svd_wrapper<-function(x=NULL, y=NULL, ind,...) {
   ## what samples: get bootstrapped data set 
    # ind contains the sample population;
    # build the bootstrap vector with repetitions
    wh<-c();
    for(i in 1:max(ind)) wh<-c(wh,which(ind>=i))

   ## what features: remove constant variables
    wp<-which(unlist(lapply(1:ncol(x),function(A) length(unique(x[wh,A]))>1)))

    numLambdas=11;
    ridgeValues=c(-5:5);
    lambdas=10^(ridgeValues); #vector of possible lambdas
    
    bestRidgeValue=ridgeValues[1];

    # the current sample and feature subset for the node
    X<-scale(x[wh,wp], scale=F)

    lbd<-rep(0, numLambdas); #vector for number of matches per lambda
    lab<-c(-1,1)[y+1] # convert two class case into fit to -1 and 1

   ## if there are no oobs left there is no need for testing different lambdas
    if(sum(ind==0)==0)
    {
      z<-gen.ridge(x=t(t(X)), y=lab[wh], lambda=lambdas[1]);
    }
    else #test different lambdas
    {
      #save objects
      objs<-vector("list",numLambdas);
      
      # the oobs with the current features
      Xt<-scale(x[ind==0,wp],  center=attr(X,"scaled:center"), scale=F)

      #SVD
      A=t(t(X));
      B=t(t(Xt));
#      s<-svd(x=A);
      s<-svd(x=A, nu=0, nv=ncol(A));
      #precalculate everything we can
      vH_Ab=Conj(t(s$v)) %*% Conj(t(A)) %*% lab[wh];

      #make sure we don't run into trouble because of
      #low rank
      sigmas<-rep(0,ncol(A));
      sigmas[1:length(s$d)]<-s$d;
      sigmasSquared<-sigmas^2;
      # loop over lambda
      for(l in 1:numLambdas){
	  
	  #the regularised pseudo inverse
	  svdCoeffs<-s$v%*%diag(1/((sigmasSquared+lambdas[l]) ) ) %*% vH_Ab;
	  objs[[l]]$coefficients<-svdCoeffs;
	  
	  #calc distance to hyperplane
	  p<-B %*% svdCoeffs;

	  # find optimal lambda from (local) oob: check how many classes match
	  lbd[l]<-sum((p>0)*1==y[ind==0])

	  #for comparison
	  #tmp<-gen.ridge(x=t(t(X)), y=lab[wh], lambda=lambdas[l])
	  #p<-predict(tmp, t(t(Xt)));
      }
#browser()
      # find (first) best lambda (with best matching)
      wm<-which.max(as.numeric(lbd))[1]
      # get again the stats for the best lambda
      #z<-gen.ridge(x=t(t(X)), y=lab[wh], lambda=lambdas[wm])
      z<-objs[[wm]];
      bestRidgeValue=ridgeValues[wm];
    }

    COEF<-z$coefficients[,1]

    v_norm=sqrt(COEF%*%COEF);
    z$coef<-rep(0,ncol(x))
    z$coef[wp]<-COEF/v_norm;

    z$signif<-rep(0,ncol(x))
    z$signif[wp]<-0;

    z$ridge_value = bestRidgeValue;
    z$lab<-"ridge"

    return(z);
}



##--------------------------------------------
## svm wrapper
##--------------------------------------------

svm_wrapper<-function(x=NULL, y=NULL, ind,...) {

   ## what samples: get bootstrapped data set 
   wh<-c();
   for(i in 1:max(ind)) wh<-c(wh,which(ind>=i))

   ## what features: remove constant variables
   wp<-which(unlist(lapply(1:ncol(x),function(A) length(unique(x[wh,A]))>1)))


    # loop over lambda
    lbd<-c()
     
    X<-scale(x[wh,wp])
    Xt<-scale(x[ind==0,wp], scale=attr(X,"scaled:scale"), center=attr(X,"scaled:center"))

    for(l in 1:3){
        
        obj<-try(svm(y=as.factor(y[wh]),x=as.matrix(X),kernel='linear', scale=T, tolerance=1, cost=10^(c(-3,-1,1))[l]),T)
        
        # find optimal lambda from oob
        if(sum(ind==0)>1)
            lbd[l]<-sum(predict(obj, newdata=Xt)==y[ind==0])

    }
    wm<-which.max(as.numeric(lbd))[1]
    z<-try(svm(y=as.factor(y[wh]),x=as.matrix(X),kernel='linear', scale=T, cost=10^(c(-3,-1,1))[wm],tolerance=1),TRUE)

    if(class(z)!="try-error")
    {
        COEF<-as.numeric(t(z$coefs)%*%as.matrix((X)[z$index,]))/attr(X,"scaled:scale")  

        v_norm=sqrt(COEF%*%COEF);
        z$coef<-rep(0,ncol(x))
        z$coef[wp]<-COEF/v_norm;

        z$signif<-rep(0,ncol(x))
        z$signif[wp]<-0;

        z$lab<-"svm"
        z$ridge_value = lbd
    }
    return(z);
}


##--------------------------------------------
## random coefficients
##--------------------------------------------

rnd_wrapper<-function(x=NULL, y=NULL, ind,...) {

   ## what samples: get bootstrapped data set 
   wh<-c();
   for(i in 1:max(ind)) wh<-c(wh,which(ind>=i))

   ## what features: remove constant variables
   wp<-which(unlist(lapply(1:ncol(x),function(A) length(unique(x[wh,A]))>1)))

    z<-list()

    X<-scale(x[wh,wp])

    COEF<-rnorm(ncol(x))

    v_norm=sqrt(COEF%*%COEF)/attr(X,"scaled:scale");
    z$coef<-COEF/v_norm;

    z$signif<-rep(0,ncol(x))
    z$ridge_value = 0    
    z$lab<-"rnorm"
    
    return(z);
}
