###############################################################################
## search.hillclimber.R --- 
## Author          : Fraser Lewis
## Last modified   : 03/10/2014
###############################################################################

## fit a given DAG to data using the distribution given
search.hillclimber <- function(score.cache=NULL,num.searches=1,seed=0, verbose=FALSE, timing.on=TRUE,start.dag=NULL,trace=FALSE,support.threshold=0.5,create.graph=FALSE,
                               dag.retained=NULL){


   data.df<-score.cache$data.df;## n.b. this might be adjusted from original data.df depending on the group variable
   ## check that all nodes are included in the cache - easy to forget if parallelising
   if(length(unique(score.cache$children))!=dim(data.df)[2]){stop("all nodes must be included in the cache!");}
   ## check valid arguments  
   if(!is.numeric(num.searches) || length(num.searches)!=1){stop("invalid num.searches");}
   if(seed<0){stop("seed must be a positive integer");}
   if(!is.logical(verbose)){ stop("invalid verbose - must be logical");}     
   if(!is.logical(timing.on)){ stop("invalid timing.on");}
   if(!is.logical(trace)){ stop("invalid trace value must be logical");}  
   if(!is.numeric(support.threshold) || support.threshold<0 || support.threshold>1){stop("invalid support.threshold value");}
   if(trace==TRUE && !interactive()){stop("asked for graphical trace but not running interactively");}
   ## create a random network from which to commence searching. This is a random combination chosen from all valid
   ## parent combination i.e. in the node cache and then checked for cyclicity

   if(!is.null(dag.retained)){check.valid.dag(dag.retained,data.df=score.cache$data.df,is.ban.matrix=FALSE,group.var=NULL);
      } else {dag.retained<-check.valid.dag(dag.retained,data.df=score.cache$data.df,is.ban.matrix=FALSE,group.var=NULL);} ##if null just create empty mat and return


   ## check that the cache has no invalid nodes - these should be removed
   myNA<-which(is.na(score.cache$mlik));
   if(length(myNA)>0){ if (verbose) cat("### NOTE: the score.cache has missing values in mlik - assigning to -.Machine$double.xmax ###\n");
   score.cache$mlik[myNA]<- -.Machine$double.xmax;## the most negative number possible, i.e. -infinity
   }

  use.start.dag<-FALSE;
  ## user supplied start network
  if(!is.null(start.dag)){ ## must also have data passed
                          if(is.null(data.df)){stop("must supply data.df when using explicit start.dag");}
                          if(num.searches!=1){stop("num.searches must equal 1 when using explicit start dag - otherwise an identical search is repeated\n => the final DAG identified depends only on the initial search location\n");}
                          check.valid.dag(dag.m=start.dag,data.df=data.df, is.ban.matrix=FALSE,NULL);
                          use.start.dag<-TRUE;}

  cache.defn<-as.integer(score.cache[["node.defn"]]);## into one long vector filled by col
  children<-as.integer(score.cache[["children"]]);
  nodescores<-as.double(score.cache[["mlik"]]);
  numVars<-as.integer(dim(score.cache[["node.defn"]])[2]);## number of variables in a DAG
  numRows<-as.integer(dim(score.cache[["node.defn"]])[1]);## total number of different parent combinations
  
  numparents.per.node<-as.integer(table(children));## number of parent combinations per variable

  ## NOTE - NaN's/NAs - these are dealt with in C utility.c make_nodecache
  ################################################################################################
  ### Want graphical output
  ### via Graphviz [and Cairo - commented out by Marta Pittavino on 03/10/2014]
  ################################################################################################
  if(trace){
      if(!requireNamespace("Rgraphviz", quietly = TRUE)){stop("library Rgraphviz is not available!\nRgraphviz is available as part of the bioconductor project - see http://www.bioconductor.org/install\n");}
  ##  if(!requireNamespace("Cairo", quietly = TRUE)){stop("library Cairo is not available!\nThis can be obtained from CRAN.");}
      
      num.traced.searches<-1;## fix at one - one plot output per search
      
      graphicsetup();## deals with opening or finding a Cairo window
      
      for(searchIndex in 1:num.searches){ if (verbose)
           cat("processing search...", searchIndex,"\n");
           if(searchIndex==1){## first time search   
         
           ## 1. run search to get first localloptimal DAG
           res<-.Call("searchhill",children,cache.defn,nodescores,numVars,numRows, numparents.per.node,as.integer(seed),
                          verbose,as.integer(timing.on), as.integer(use.start.dag),as.integer(start.dag),as.integer(num.traced.searches),
                          as.integer(dag.retained),
               PACKAGE="abn" ## uncomment to load as package not shlib
              )
           for(i in 2:length(res)){colnames(res[[i]])<-colnames(score.cache[[2]]);rownames(res[[i]])<-colnames(score.cache[[2]]);} ## just for names
           
           init.score<-rep(NA,num.searches);init.score[searchIndex]<-res[[1]][1];
           final.score<-rep(NA,num.searches); final.score[searchIndex]<-res[[1]][2];         
           init.dag<-list();init.dag[[searchIndex]]<-res[2][[1]]; 
           final.dag<-list();final.dag[[searchIndex]]<-res[3][[1]]; 

           con.dag<-res[2+1][[1]];## the best DAG from a single search
           mygraph<-new("graphAM",adjMat=t(con.dag),edgemode="directed");
           mygraph.prev<-new("graphAM",adjMat=t(con.dag),edgemode="directed");
           if(num.searches>1){
           plot(mygraph,main=paste("Consensus DAG - at ",100*support.threshold,"% support - still searching....",sep=""));
           }
           } else {
                   ## have the second or more iteration and so compute the consensus
                   seed<-seed+1;## increment seed
                   ## 1. run search to get first localloptimal DAG
                   res<-.Call("searchhill",children,cache.defn,nodescores,numVars,numRows, numparents.per.node,as.integer(seed),
                                           verbose,as.integer(timing.on), as.integer(use.start.dag),as.integer(start.dag),as.integer(num.traced.searches),
                                           as.integer(dag.retained),
                              PACKAGE="abn" ## uncomment to load as package not shlib
                              )
                         for(i in 2:length(res)){colnames(res[[i]])<-colnames(score.cache[[2]]);rownames(res[[i]])<-colnames(score.cache[[2]]);}

                   init.score[searchIndex]<-res[[1]][1];
                   final.score[searchIndex]<-res[[1]][2];
                   init.dag[[searchIndex]]<-res[2][[1]]; 
                   final.dag[[searchIndex]]<-res[3][[1]]; 

           con.dag<-con.dag+res[2+1][[1]];## the best DAG from a single search and then add to existing to get freq counts
           #cur.score<-res[[1]][2];
           con.dag.binary<-ifelse(con.dag>=ceiling(searchIndex*support.threshold),1,0);## drop all arcs that are not in a majority >50% of searches so far run
           mygraph<-initialize(mygraph,adjMat=t(con.dag.binary),edgemode="directed");##update graph
           if(!identical(mygraph,mygraph.prev)){## if not identical then re-do plot
                                              plot(mygraph,main=paste("Consensus DAG - at ",100*support.threshold,"% support - still searching....",sep=""));
                                              mygraph.prev<-initialize(mygraph,adjMat=t(con.dag.binary),edgemode="directed");## update
           } else{ ## do nothing but same no change
                   if (verbose) cat("no change in consensus\n");}
                      } ## else for searchIndex==1 
           }
           
          ## do a final plot to avoid refresh issue
          plot(mygraph,main=paste("Consensus DAG - at ",100*support.threshold,"% support after ",num.searches, " searches",sep=""));
          ## return results in similar format to batch method below
       if(num.searches==1){con.dag.binary<-con.dag;}

      if(create.graph){## create graph object part 

      mygraph<-new("graphAM",adjMat=t(con.dag.binary),edgemode="directed");

      return(list(init.score=init.score,final.score=final.score,init.dag=init.dag,final.dag=final.dag,consensus=con.dag.binary,support.threshold=support.threshold,graph=mygraph));
                
      } else { 
              return(list(init.score=init.score,final.score=final.score,init.dag=init.dag,final.dag=final.dag,consensus=con.dag.binary,support.threshold=support.threshold));     
             }
          


  } else { 
  ################################################################################################
  ### No graphical output - all inside C
  ### 
  ################################################################################################
  ## get a random DAG from which to start - N.B. last argument is verbose
  res<-.Call("searchhill",children,cache.defn,nodescores,numVars,numRows, numparents.per.node,as.integer(seed),
                          verbose,as.integer(timing.on), as.integer(use.start.dag),as.integer(start.dag),as.integer(num.searches),
                          as.integer(dag.retained),
              PACKAGE="abn" ## uncomment to load as package not shlib
              ) 

    ## results format is res[[1]] is a vector of network scores in order: init score1, final score1, init score2, final score2,....etc
    ## with res[[2]] = init network 1, res[[3]] final network 1, res[[4]] init network2, res[[5]] init network2,... etc
    for(i in 2:length(res)){colnames(res[[i]])<-colnames(score.cache[[2]]);rownames(res[[i]])<-colnames(score.cache[[2]]);}
    ## now re-organise res into something easier to analyse - a list of three lists
    ## list 1 - vector of scores and list of initial matrices, list 2 -vector of score and list of final matrices
    scores<-res[[1]];
     init.indexes<-seq(1,2*num.searches,by=2);
     fin.indexes<-seq(2,2*num.searches,by=2);
                     
    init.scores<-scores[init.indexes];## scores from initial networks
    fin.scores<-scores[fin.indexes];
    init.mat<-res[init.indexes+1];##offset for score vector
    fin.mat<-res[fin.indexes+1];##offset for score vector
    rm(res);

    con.dag<-fin.mat[[1]];
    if(num.searches>1){for(i in 2:num.searches){con.dag<-con.dag+fin.mat[[i]];}} ## get total arc freq
    ## now make binary according to threshold passed
    con.dag.binary<-ifelse(con.dag>=ceiling(num.searches*support.threshold),1,0);
    #########################################################
    ## create graph object part
    #########################################################
      if(create.graph){
      if(!requireNamespace("Rgraphviz", quietly = TRUE)){stop("library Rgraphviz is not available!\nRgraphviz is available as part of the bioconductor project - see http://www.bioconductor.org/install\nRgraphviz is required is create.graph=TRUE");}
      mygraph<-new("graphAM",adjMat=t(con.dag.binary),edgemode="directed");
      return(list(init.score=init.scores,final.score=fin.scores,init.dag=init.mat,final.dag=fin.mat,consensus=con.dag.binary,support.threshold=support.threshold,graph=mygraph));
     } else {
             return(list(init.score=init.scores,final.score=fin.scores,init.dag=init.mat,final.dag=fin.mat,consensus=con.dag.binary,support.threshold=support.threshold));
               }
 
  } ## end of C part

}
