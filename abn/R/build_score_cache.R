###############################################################################
## buildscorecache.R --- 
## Author          : Fraser Lewis
## Last modified   : 02/12/2012
###############################################################################

## fit a given DAG to data
buildscorecache <- function(data.df=NULL, data.dists=NULL, group.var=NULL,cor.vars=NULL,
                            dag.banned=NULL, dag.retained=NULL,
                            max.parents=NULL,
                            which.nodes=NULL,defn.res=NULL,dry.run=FALSE,
                            max.mode.error=10,
                            verbose=FALSE,centre=TRUE,mean=0, prec=0.001,loggam.shape=1,
                            loggam.inv.scale=5e-05, max.iters=100,
                            epsabs=1e-7,error.verbose=FALSE,
                            epsabs.inner=1e-6,max.iters.inner=100,
                            finite.step.size=1e-7,
                            hessian.params=c(1E-04,1E-02),max.iters.hessian=10,
                            max.hessian.error=5E-01,factor.brent=1E+02, maxiters.hessian.brent=100,num.intervals.brent=100
                  ){

   ## simple check
   if(!is.null(cor.vars) && is.null(group.var)){stop("must specify the group variable!");}  
   
   ## another simple check   
   if(max.mode.error>100 || max.mode.error<0){stop("max.mode.error is a % and must be [0,100]!");}
   
   force.method<-"notset";
   if(max.mode.error==0){force.method<-"C";}
   if(max.mode.error==100){force.method<-"INLA";}

   ## run series of checks on the arguments
   mylist<-check.valid.data(data.df,data.dists,group.var);## return a list with entries bin, gaus, pois, ntrials and exposure

   ## run a series of common sense checks on the banned DAG and retain DAG
      if(!is.null(dag.banned)){check.valid.dag(dag.banned,data.df=data.df,is.ban.matrix=TRUE,group.var=group.var);
      } else {dag.banned<-check.valid.dag(dag.banned,data.df=data.df,is.ban.matrix=TRUE,group.var=group.var);} ##if null just create empty mat and return

      if(!is.null(dag.retained)){check.valid.dag(dag.retained,data.df=data.df,is.ban.matrix=FALSE,group.var=group.var);
      } else {dag.retained<-check.valid.dag(dag.retained,data.df=data.df,is.ban.matrix=FALSE,group.var=group.var);} ##if null just create empty mat and return

    ## check max.parents is a list with suitable entries      
      if(is.null(defn.res)){## not supplying custom parent sets
        max.parents<-check.valid.parents(data.df,max.parents,group.var);## returns an integer vector of the same length of number of nodes
      } else {## providing custom parent set and so max.parents is irrelevant
              max.parents<-max(apply(defn.res[["node.defn"]],1,sum));}## note this is just a dummy and not actually used anywhere!
      ## check retain does not ask for more arcs to be retained than allowed in max.parents
      max.retain<-apply(dag.retained,1,sum);## number of parents per node to retain
      if(length(which( (max.retain>max.parents) == TRUE))>0){stop("dag.retained is inconsistent with max.parents!");}

      ## check that arcs than are banned are also not retained      
      if(length(which(which(as.integer(dag.banned)==1)%in%which(as.integer(dag.retained)==1)==TRUE))>0){stop("dag.banned and dag.retained are inconsistent!");}
      
      ## check which.nodes is sensible
      if(is.null(defn.res)){which.nodes<-check.which.valid.nodes(data.df,which.nodes,group.var);
      } else {## have user supplied children and parent combinations 
              which.nodes<-unique(defn.res$children);}

   ## check grouping variables
   list.group.var<-check.valid.groups(group.var,data.df,cor.vars);## returns ammended data.df and suitable variables
   grouped.vars<-list.group.var$grouped.vars;## int vect of variables to be treated as grouped indexed from 1
   group.ids<-list.group.var$group.ids;## int vector of group membership ids
   data.df<-list.group.var$data.df;## this has removed the grouping variable from data.df

   ## check split.vars - gives the variables which are multinomial but split into binary  
   ## do this check after removing the grouping variable from data.df (if that was necessary)
    #  if(!is.null(split.vars)){if(length(split.vars)!=dim(data.df)[2]){stop("split.vars is the wrong length!");}
    #                           if(max(split.vars)>=dim(data.df)[2]){stop("split.vars seems in wrong format - should be integers from 1 through to number of unique variables");}
    #                           }

   ## coerce binary factors to become 0/1 integers - the 0 is based on the first entry in levels()
   if(!is.null(mylist$bin)){## have at least one binary variable
      for(i in mylist$bin){data.df[,i]<-as.numeric(data.df[,i])-1;}
      }

   ## standardize gaussian variables to zero mean and sd=1
   if(centre && !is.null(mylist$gaus)){## have at least one gaussian variable
        for(i in mylist$gaus){data.df[,i]<-(data.df[,i]-mean(data.df[,i]))/sd(data.df[,i]);}
      }

   ## coerce all data points to doubles
   for(i in 1:dim(data.df)[2]){data.df[,i]<-as.double(data.df[,i]);}## coerce ALL cols to double equivalents

   ## check for type of grouped variable
   # if(!is.null(group.var)){
   #   for(i in grouped.vars){## for each variable to be treated as grouped 
   #     if(data.dists[[i]]!="binomial"){## 
   #                        stop("currently grouped variables must be binary only");}
   #   }}     

  ##
  var.types<-get.var.types(data.dists); ## get distributions in terms of a numeric code

 
########################################################################################
## All checking over
## get to here we have suitable data/variables and now generate all the parent combinations
########################################################################################
 ## down to here we have all the data correct and now call C buildnodecache() to create all the node definitions. 
     if(is.null(defn.res)){
      ## pass to C the number (number_of_nodes,banned_arc_as_vector, retain_arcs_as_vector, max_parents_as_vector
      res<-.Call("buildcachematrix",dim(dag.banned)[1],as.integer(dag.banned),as.integer(dag.retained), max.parents, which.nodes
              ,PACKAGE="abn" ## uncomment to load as package not shlib
              ); 

     defn.res<-list();
     defn.res[["children"]]<-res[[1]];
     defn.res[["node.defn"]]<-matrix(data=res[[2]],byrow=TRUE,ncol=dim(data.df)[2]);
               colnames(defn.res[["node.defn"]])<-names(data.df);
      rm(res);
      
    # if(!is.null(split.vars)){adjust.for.split.variables(split.vars,defn.res);}

      } else { ## some check since user has supplied defn.res
               if(!is.list(defn.res)){stop("defn.res must be a list");}
               if(length(defn.res)!=2){stop("defn.res must have two entries");}
               if(!(is.vector(defn.res[[1]]) && is.matrix(defn.res[[2]]))){stop("defn.res is wrong format");}
               if(!(max(defn.res[[2]])==1 && min(defn.res[[2]])==0)){stop("defn.res is wrong format - must only be 0,1 in node definitions");} 
              }
      
      if(dry.run){## don't do any computation just return the node definitions
                  if (verbose) cat("No computation - returning only the node combinations\n"); 
                  return(defn.res);}

### loop through each node and find out what kind of model is to be fitted and then pass to appropriate
### separate call for each individual node-parentcombination
 total.rows<-dim(defn.res[["node.defn"]])[1];
 defn.res[["mlik"]]<-rep(NA,total.rows);
 defn.res[["error.code"]]<-rep(NA,total.rows);#res.list[["error.code"]]<-NULL;
 defn.res[["hessian.accuracy"]]<-rep(NA,total.rows);#NULL;res.list[["hessian.accuracy"]]<-NULL;
 defn.res[["used.INLA"]]<-rep(NA,total.rows);
 dag.m<-matrix(rep(NA,dim(data.df)[2]^2),ncol=dim(data.df)[2]);## just a dummy so code easier matches fitabn()
 colnames(dag.m)<-rownames(dag.m)<-names(data.df);
###########################################################
## Iterate over each node in the DAG separately
###########################################################
#cat("########################################################\n");
#cat("###### Fitting cache to data\n");
row.num<-1;
for(child in defn.res[["children"]]){## for each child in the cache
FAILED<-FALSE;## to catch any crashes...
###########################################################
########################################################### 
#cat("###### Processing...",row.num," of ",total.rows,"\n");
                           dag.m[,]<-0;## reset to empty
                           dag.m[child,]<-defn.res[["node.defn"]][row.num,];## set parent combination                          
                           orig.force.method<-NULL;
                           used.inla<-TRUE;
                           ####################################################
                           ### First case is the node a GLM
                           ####################################################
                           if( !(child%in%grouped.vars)){## only compute option here is C since fast and INLA slower and less reliable
                               if(force.method=="notset" || force.method=="C"){
                                 #cat("Using internal code (Laplace glm)\n");
                                 r<-try(res.c <- .Call("fit_single_node",
                                                data.df,
                                                as.integer(child),## childnode
                                                as.integer(dag.m[child,]),## parent combination
                                                as.integer(dim(dag.m)[1]),## number of nodes/variables
                                                as.integer(var.types),## type of densities
                                                as.integer(sum(dag.m[child,])),## max.parents
                                                as.double(mean),as.double(1/sqrt(prec)),as.double(loggam.shape),as.double(1/loggam.inv.scale),
                                                as.integer(max.iters),as.double(epsabs),
                                                as.integer(verbose),as.integer(error.verbose),
                                                as.integer(grouped.vars-1),## int.vector of variables which are mixed model nodes -1 for C
                                                group.ids,## group memberships - note indexed from 1
                                                as.double(epsabs.inner),
                                                as.integer(max.iters.inner),
                                                as.double(finite.step.size),
                                                as.double(hessian.params),
                                                as.integer(max.iters.hessian),
                                                as.integer(0), ## modes only - false here as only applies to glmms
                                                as.double(max.hessian.error),## Not applicable 
                                                as.double(factor.brent),## Not applicable 
                                                as.integer(maxiters.hessian.brent),## Not applicable 
                                                as.double(num.intervals.brent)## Not applicable 
                                  ,PACKAGE="abn" ## uncomment to load as package not shlib
                                                 ));#print(res);
                          if(length(attr(r,"class")>0) && attr(r,"class")=="try-error"){ if (verbose) cat("## !!! Laplace approximation failed\n");
                                                                           FAILED<-TRUE;} 
                           used.inla<-FALSE;## flip
                           } else {## use INLA for glm
                                  if(!requireNamespace("INLA", quietly = TRUE)){stop("library INLA is not available!\nR-INLA is available from http://www.r-inla.org/download\nAfter installation please use inla.upgrade() to get the latest version (this is required)");}
                                  mean.intercept<-mean;## use same as for rest of linear terms 
                                  prec.intercept<-prec;## use same as for rest of linear terms
                                  #cat("Using INLA (glm)\n");
                                         res.inla<-calc.node.inla.glm(child,
                                                                      dag.m,
                                                                      data.df,
                                                                      data.dists,
                                                                      rep(1,dim(data.df)[1]),## ntrials
                                                                      rep(1,dim(data.df)[1]),## exposure  
                                                                      TRUE, mean.intercept, prec.intercept, mean, prec,loggam.shape,loggam.inv.scale,verbose);
                              
                                  if(is.logical(res.inla)){ if (verbose) cat("INLA failed....so reverting to internal code\n");
                                                           r<-try(res.c <- .Call("fit_single_node",
                                                                          data.df,
                                                                          as.integer(child),## childnode
                                                                          as.integer(dag.m[child,]),## parent combination
                                                                          as.integer(dim(dag.m)[1]),## number of nodes/variables
                                                                          as.integer(var.types),## type of densities
                                                                          as.integer(sum(dag.m[child,])),## max.parents
                                                                          as.double(mean),as.double(1/sqrt(prec)),as.double(loggam.shape),as.double(1/loggam.inv.scale),
                                                                          as.integer(max.iters),as.double(epsabs),
                                                                          as.integer(verbose),as.integer(error.verbose),
                                                                          as.integer(grouped.vars-1),## int.vector of variables which are mixed model nodes -1 for C
                                                                          group.ids,## group memberships - note indexed from 1
                                                                          as.double(epsabs.inner),
                                                                          as.integer(max.iters.inner),
                                                                          as.double(finite.step.size),
                                                                          as.double(hessian.params),
                                                                          as.integer(max.iters.hessian),
                                                                          as.integer(0), ## modes only - false here as only applies to glmms
                                                                          as.double(max.hessian.error),## Not applicable 
                                                                          as.double(factor.brent),## Not applicable 
                                                                          as.integer(maxiters.hessian.brent),## Not applicable 
                                                                          as.double(num.intervals.brent)## Not applicable 
                                                                          ,PACKAGE="abn" ## uncomment to load as package not shlib
                                                                         ));#print(res);
                                                         if(length(attr(r,"class")>0) && attr(r,"class")=="try-error"){ if (verbose) cat("## !!! Laplace approximation failed\n");
                                                                           FAILED<-TRUE;}  
                                                             used.inla<-FALSE;## flip        
                                            } ## INLA failed                                        
                           
                           } ## use INLA
                           ###########################################################
                           ###########################################################
                           ## End of GLM node
                           ########################################################### 
                         
                           } else {
                           ###########################################################
                           ## Have a GLMM node
                           ########################################################### 
                                  ## have a glmm, so two options, INLA or C
                                  if(force.method=="notset" || force.method=="INLA"){##  
                                  if(!requireNamespace("INLA", quietly = TRUE)){stop("library INLA is not available!\nR-INLA is available from http://www.r-inla.org/download\nAfter installation please use inla.upgrade() to get the latest version (this is required)");}
                                  mean.intercept<-mean;## use same as for rest of linear terms 
                                  prec.intercept<-prec;## use same as for rest of linear terms
                                  res.inla<-calc.node.inla.glmm(child,
                                                            dag.m,
                                                            data.frame(data.df,group=group.ids),
                                                            data.dists,
                                                            rep(1,dim(data.df)[1]),## ntrials
                                                            rep(1,dim(data.df)[1]),## exposure 
                                                            TRUE,## always compute marginals - since only way to check results
                                                            mean.intercept, prec.intercept, mean, prec,loggam.shape,loggam.inv.scale,verbose); 
                                   #return(res);

                                   ## CHECK FOR INLA CRASH
                                   if(is.logical(res.inla)){ if (verbose) cat("INLA failed....so reverting to internal code\n");
                                                            orig.force.method<-force.method;## save original
                                                            force.method="C"; ## Easiest way is just to force C for this node
                                   } else {
                                   res.inla.modes<-getModeVector(list.fixed=res.inla$marginals.fixed,list.hyper=res.inla$marginals.hyperpar);}
                                   #cat("check INLA modes against C at node ",rownames(dag.m)[child],"\n");
                                   
                                   }
                                   #cat("fit a glmm at node ",rownames(dag.m)[child],"using C\n");
                               if(force.method=="notset"){
                                  r<-try(res.c <- .Call("fit_single_node",
                                            data.df,
                                            as.integer(child),## childnode
                                            as.integer(dag.m[child,]),## parent combination
                                            as.integer(dim(dag.m)[1]),## number of nodes/variables
                                            as.integer(var.types),## type of densities
                                            as.integer(sum(dag.m[child,])),## max.parents
                                            as.double(mean),as.double(1/sqrt(prec)),as.double(loggam.shape),as.double(1/loggam.inv.scale),
                                            as.integer(max.iters),as.double(epsabs),
                                            as.integer(verbose),as.integer(error.verbose),
                                            as.integer(grouped.vars-1),## int.vector of variables which are mixed model nodes -1 for C
                                            group.ids,## group memberships - note indexed from 1
                                            as.double(epsabs.inner),
                                            as.integer(max.iters.inner),
                                            as.double(finite.step.size),
                                            as.double(hessian.params),
                                            as.integer(max.iters.hessian),
                                            as.integer(1), ## turn on ModesONLY
                                            as.double(max.hessian.error),## Not applicable 
                                            as.double(factor.brent),## Not applicable 
                                            as.integer(maxiters.hessian.brent),## Not applicable 
                                            as.double(num.intervals.brent)## Not applicable 
                                  ,PACKAGE="abn" ## uncomment to load as package not shlib
                                                 ));#print(res);  
                            
                             if(length(attr(r,"class")>0) && attr(r,"class")=="try-error"){ if (verbose) cat("## !!! Laplace approximation failed\n");
                                                                           FAILED<-TRUE;}     
                            if(!FAILED){              
                            res.c.modes<-res.c[[1]][-c(1:3)];## remove mlik - this is first entry, and error code and hessian accuracy
                            res.c.modes<-res.c.modes[which(res.c.modes!=.Machine$double.xmax)];## this discards all "empty" parameters 
                            ## get difference in modes proportion relative to C
                            diff.in.modes<-(res.inla.modes-res.c.modes)/res.c.modes;
                            error.modes<-max(abs(diff.in.modes));
                            }
                            } ## end of notset

                           if( !FAILED && force.method=="C" || (force.method=="notset" && error.modes>(max.mode.error/100))){ ## INLA might be unreliable so use C (slower)
                           
                               if(force.method=="notset"){ if (verbose) cat("Using internal code (Laplace glmm)\n=>max. abs. difference (in %) with INLA is ");
                               if (verbose) cat(formatC(100*error.modes,format="f",digits=1)," and exceeds tolerance\n");} else { if (verbose) cat("Using internal code (Laplace glmm)\n");}

                               r<-try(res.c <- .Call("fit_single_node",
                                            data.df,
                                            as.integer(child),## childnode
                                            as.integer(dag.m[child,]),## parent combination
                                            as.integer(dim(dag.m)[1]),## number of nodes/variables
                                            as.integer(var.types),## type of densities
                                            as.integer(sum(dag.m[child,])),## max.parents
                                            as.double(mean),as.double(1/sqrt(prec)),as.double(loggam.shape),as.double(1/loggam.inv.scale),
                                            as.integer(max.iters),as.double(epsabs),
                                            as.integer(verbose),as.integer(error.verbose),
                                            as.integer(grouped.vars-1),## int.vector of variables which are mixed model nodes -1 for C
                                            group.ids,## group memberships - note indexed from 1
                                            as.double(epsabs.inner),
                                            as.integer(max.iters.inner),
                                            as.double(finite.step.size),
                                            as.double(hessian.params),
                                            as.integer(max.iters.hessian),
                                            as.integer(0), ## turn on ModesONLY
                                            as.double(max.hessian.error),## Not applicable 
                                            as.double(factor.brent),## Not applicable 
                                            as.integer(maxiters.hessian.brent),## Not applicable 
                                            as.double(num.intervals.brent)## Not applicable 
                                  ,PACKAGE="abn" ## uncomment to load as package not shlib
                                                 ));
                                  if(length(attr(r,"class")>0) && attr(r,"class")=="try-error"){ if (verbose) cat("## !!! Laplace approximation failed\n");
                                                                           FAILED<-TRUE;} 
                                   
                                   used.inla<-FALSE;## flip
                                   
                                 } else { if (verbose) cat("Using INLA (glmm)\n");
                                         }## end of if inla bad
                            } ## end of if GLMM
                            ###########################################################
                            ## End of GLMM node
                            ########################################################### 
                           
                            ###########################################################
                            ## End of all external computations
                            ###########################################################
                            ## computation for current node is all done so sort out the 
                            ## output into nicer form and give labels
                            ###########################################################
                            child.name<-colnames(dag.m)[child];
  
                            if(!FAILED){                            
                                if(used.inla==FALSE){## organize output from C
                                defn.res[["mlik"]][row.num]<-res.c[[1]][1];            
                                defn.res[["error.code"]][row.num]<-res.c[[1]][2];
                                defn.res[["hessian.accuracy"]][row.num]<-res.c[[1]][3];
                                defn.res[["used.INLA"]][row.num]<-FALSE;
                                } else {
                                       ## organize output from INLA
                                       defn.res[["mlik"]][row.num]<-res.inla$mlik[2];## [2] is for Gaussian rather than Integrated estimate 
                                       defn.res[["error.code"]][row.num]<-NA;## not available from INLA
                                       defn.res[["hessian.accuracy"]][row.num]<-NA;## not available from INLA
                                       defn.res[["used.INLA"]][row.num]<-TRUE;
                                        } 
                            } else {## FAILED
                                   defn.res[["mlik"]][row.num]<-NA;            
                                   defn.res[["error.code"]][row.num]<-2;## model could not be fitted
                                   defn.res[["hessian.accuracy"]][row.num]<-NA;
                                   defn.res[["used.INLA"]][row.num]<-FALSE;
                                      }
                              
                            if(!is.null(orig.force.method)){force.method<-orig.force.method;} ## reset force.method after INLA crash
                            ############################################################
                            ## Finished with current node
                            ############################################################

row.num<-row.num+1;
 } ## end of nodes loop
                           
      #################################################
      ## Further tidy up of results across all nodes
      #################################################
      
      defn.res[["error.code.desc"]]<-as.character(defn.res[["error.code"]]);
      defn.res[["error.code.desc"]]<-ifelse(defn.res[["error.code.desc"]]==0,"success",defn.res[["error.code"]]);
      defn.res[["error.code.desc"]]<-ifelse(defn.res[["error.code.desc"]]==1,"warning: mode results may be unreliable (optimiser terminated unusually)",defn.res[["error.code.desc"]]);
      defn.res[["error.code.desc"]]<-ifelse(defn.res[["error.code.desc"]]==2,"error - logscore is NA - model could not be fitted",defn.res[["error.code.desc"]]);
      defn.res[["error.code.desc"]]<-ifelse(defn.res[["error.code.desc"]]==4,"warning: fin.diff hessian estimation terminated unusually ",defn.res[["error.code.desc"]]);

      defn.res[["data.df"]]<-data.df;
      defn.res[["group.ids"]]<-group.ids;

#cat("########End of cache building ##########################\n");
#cat("########################################################\n");
return(defn.res);

}

