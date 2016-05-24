###############################################################################
## fitabn.R --- 
## Author          : Fraser Lewis
## Last modified   : 02/12/2012
## Last modified   : 29/09/2014 by Marta Pittavino, just renamed.
###############################################################################

## fit a given DAG to data
fitabn <- function(dag.m=NULL, data.df=NULL, data.dists=NULL, group.var=NULL,cor.vars=NULL,create.graph=FALSE,compute.fixed=FALSE,
                   mean=0, prec=0.001,
                   loggam.shape=1,loggam.inv.scale=5e-05,verbose=FALSE, centre=TRUE,
                   max.mode.error=10,
                   max.iters=100,epsabs=1e-7,error.verbose=FALSE,epsabs.inner=1e-6,max.iters.inner=100,
                   finite.step.size=1E-07,hessian.params=c(1E-04,1E-02),max.iters.hessian=10,
                   max.hessian.error=1E-04,factor.brent=1E+02, maxiters.hessian.brent=10,num.intervals.brent=100,
                   min.pdf=1E-03,n.grid=100,std.area=TRUE, marginal.quantiles=c(0.025,0.25,0.5,0.75,0.975),max.grid.iter=1000,
                   marginal.node=NULL, marginal.param=NULL,variate.vec=NULL
                  ){
   
   ## simple check
   if(!is.null(cor.vars) && is.null(group.var)){stop("must specify the group variable!");}  
   
   ## another simple check   
   if(max.mode.error>100 || max.mode.error<0){stop("max.mode.error is a % and must be [0,100]!");}   
   #if(!is.null(marginal.node)){std.area<-FALSE;n.grid<-NULL;marginal.quantiles<-NULL;}
   if(is.null(max.grid.iter)){stop("max.grid.iter cannot be NULL!");}

   if(!is.null(variate.vec)){max.mode.error<-0;std.area<-FALSE;n.grid<-NULL;} ## if user supplied grid then must use C

   force.method<-"notset";
   if(max.mode.error==0){force.method<-"C";}
   if(max.mode.error==100){force.method<-"INLA";}

   ## run series of checks on the arguments
   mylist<-check.valid.data(data.df,data.dists,group.var);## return a list with entries bin, gaus, pois, ntrials and exposure

   ## run a series of checks on the DAG passed
   check.valid.dag(dag.m=dag.m,data.df=data.df,is.ban.matrix=FALSE,group.var=group.var);

   ## check grouping variables
   list.group.var<-check.valid.groups(group.var,data.df,cor.vars);## returns ammended data.df and suitable variables
   grouped.vars<-list.group.var$grouped.vars;## int vect of variables to be treated as grouped indexed from 1
   group.ids<-list.group.var$group.ids;## int vector of group membership ids
   data.df<-list.group.var$data.df;## this has removed the grouping variable from data.df

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
  #  if(!is.null(group.var)){
  #    for(i in grouped.vars){## for each variable to be treated as grouped 
  #      if(data.dists[[i]]!="binomial"){## 
  #                         stop("currently grouped variables must be binary only");}
  #    }}     

  ##
  var.types<-get.var.types(data.dists); ## get distributions in terms of a numeric code

########################################################################################
## All checking over
## get to here we have suitable data/variables and not do the actual fitting
########################################################################################

### loop through each node and find out what kind of model is to be fitted and then pass to appropriate
### separate call for each individual node
## setup some storage for entire DAG
 res.list<-list();
 error.code<-rep(NA,dim(dag.m)[1]);#res.list[["error.code"]]<-NULL;
 hessian.accuracy<-rep(NA,dim(dag.m)[1]);#NULL;res.list[["hessian.accuracy"]]<-NULL;
 mymodes<-list();
 mymargs<-list();
 INLA.marginals<-NULL;
###########################################################
## Iterate over each node in the DAG separately
###########################################################
if (verbose) cat("########################################################\n");
if (verbose) cat("###### Fitting DAG to data\n");
for(child in 1:dim(dag.m)[1]){## for each node in the DAG

###########################################################
########################################################### 
if (verbose) cat("###### Processing...Node ",rownames(dag.m)[child],"\n");                          
                           orig.force.method<-NULL;
                           used.inla<-TRUE;
                           ####################################################
                           ### First case is the node a GLM
                           ####################################################
                           if( !(child%in%grouped.vars)){## only compute option here is C since fast and INLA slower and less reliable
                           
                                if(force.method=="notset" || force.method=="C"){
                                 if (verbose) cat("Using internal code (Laplace glm)\n");
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
                                                as.integer(0),## Not applicable 
                                                as.double(max.hessian.error),## Not applicable 
                                                as.double(factor.brent),## Not applicable 
                                                as.integer(maxiters.hessian.brent),## Not applicable 
                                                as.double(num.intervals.brent)## Not applicable 
                                  ,PACKAGE="abn" ## uncomment to load as package not shlib
                                                 ));#print(res.c);
                                    if(length(attr(r,"class")>0) && attr(r,"class")=="try-error"){cat("## !!! Laplace approximation failed at node ",
                                                        rownames(dag.m)[child],
                                                        "\n## The additive formulation at this node is perhaps over-parameterized?\n",
                                                        "## Fitting the glm at this node using glm() may provide more information\n",
                                                        "## If glm() can fit this model then please send a bug report to fraseriain.lewis@uzh.ch\n",sep="");
                                                         stop("");}
                           used.inla<-FALSE;## flip
                           } else {## use INLA for glm
                                  if(!requireNamespace("INLA", quietly = TRUE)){stop("library INLA is not available!\nR-INLA is available from http://www.r-inla.org/download\nAfter installation please use inla.upgrade() to get the latest version (this is required)");}
                                  mean.intercept<-mean;## use same as for rest of linear terms 
                                  prec.intercept<-prec;## use same as for rest of linear terms 
                                  if (verbose) cat("Using INLA (glm)\n");
                                         res.inla<-calc.node.inla.glm(child,
                                                                      dag.m,
                                                                      data.df,
                                                                      data.dists,
                                                                      rep(1,dim(data.df)[1]),## ntrials
                                                                      rep(1,dim(data.df)[1]),## exposure  
                                                                      TRUE, mean.intercept, prec.intercept, mean, prec,loggam.shape,loggam.inv.scale,verbose);
                                                                     

                                          if(is.logical(res.inla)){if (verbose) cat("INLA failed....so reverting to internal code\n");
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
                                                                          as.integer(0),## Not applicable 
                                                                          as.double(max.hessian.error),## Not applicable 
                                                                          as.double(factor.brent),## Not applicable 
                                                                          as.integer(maxiters.hessian.brent),## Not applicable 
                                                                          as.double(num.intervals.brent)## Not applicable 
                                                                          ,PACKAGE="abn" ## uncomment to load as package not shlib
                                                                         ));#print(res);  
                                                        if(length(attr(r,"class")>0) && attr(r,"class")=="try-error"){cat("## !!! Laplace approximation failed at node ",
                                                        rownames(dag.m)[child],
                                                        "\n## The additive formulation at this node is perhaps over-parameterized?\n",
                                                        "## Fitting the glm at this node using glm() may provide more information\n",
                                                        "## If glm() can fit this model then please send a bug report to fraseriain.lewis@uzh.ch\n",sep="");
                                                         stop("");}

                                                             used.inla<-FALSE;## flip        
                                            } ## INLA failed

                           } ## use INLA
                           ###########################################################
                           ## End of GLM node
                           ########################################################### 
                         
                           } else {
                           ###########################################################
                           ## Have a GLMM node
                           ########################################################### 
                                  ## have a glmm, so two options, INLA or C
                                 
                                  if(force.method=="notset" || force.method=="INLA"){##  
                                  if(requireNamespace("INLA", quietly = TRUE)){stop("library INLA is not available!\nR-INLA is available from http://www.r-inla.org/download\nAfter installation please use inla.upgrade() to get the latest version (this is required)");}
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
                                   #return(res.inla); ## to return RAW INLA object
                                   ## CHECK FOR INLA CRASH
                                   if(is.logical(res.inla)){if (verbose) cat("INLA failed....so reverting to internal code\n");
                                                            orig.force.method<-force.method;## save original
                                                            force.method="C"; ## Easiest way is just to force C for this node
                                   } else {
                                   res.inla.modes<-getModeVector(list.fixed=res.inla$marginals.fixed,list.hyper=res.inla$marginals.hyperpar);}
                                   #print(res.inla.modes);stop("");
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
                                            as.integer(1),## turn on ModesONLY
                                            as.double(max.hessian.error),
                                            as.double(factor.brent),
                                            as.integer(maxiters.hessian.brent),
                                            as.double(num.intervals.brent)
                                  ,PACKAGE="abn" ## uncomment to load as package not shlib
                                                 ));#print(res);  
                              if(length(attr(r,"class")>0) && attr(r,"class")=="try-error"){cat("## !!! Laplace approximation failed at node ",
                                                        rownames(dag.m)[child],
                                                        "\n## The additive formulation at this node is perhaps over-parameterized?\n",
                                                        "## Fitting the glmm at this node using glmer() in lme4 may provide more information\n",
                                                        "## If glmer() can fit this model then please send a bug report to fraseriain.lewis@uzh.ch\n",sep="");
                                                         stop("");}
                           
                            res.c.modes<-res.c[[1]][-c(1:3)];## remove mlik - this is first entry, and error code and hessian accuracy
                            res.c.modes<-res.c.modes[which(res.c.modes!=.Machine$double.xmax)];## this discards all "empty" parameters 
                            ## get difference in modes proportion relative to C
                            diff.in.modes<-(res.inla.modes-res.c.modes)/res.c.modes;
                            error.modes<-max(abs(diff.in.modes));
                            }

                            if( force.method=="C" || (force.method=="notset" && error.modes>(max.mode.error/100))){ ## INLA might be unreliable to use C (slower)
                               
                                if(force.method=="notset"){if (verbose) cat("Using internal code (Laplace glmm)\n=>max. abs. difference (in %) with INLA is ");
                                if (verbose) cat(formatC(100*error.modes,format="f",digits=1)," and exceeds tolerance\n");} else {if (verbose) cat("Using internal code (Laplace glmm)\n");}

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
                                            as.integer(0),## turn off ModesONLY
                                            as.double(max.hessian.error),
                                            as.double(factor.brent),
                                            as.integer(maxiters.hessian.brent),
                                            as.double(num.intervals.brent) 
                                  ,PACKAGE="abn" ## uncomment to load as package not shlib
                                                 )#print(res.c);
                                   );
                                   if(length(attr(r,"class")>0) && attr(r,"class")=="try-error"){cat("## !!! Laplace approximation failed at node ",
                                                        rownames(dag.m)[child],
                                                        "\n## The additive formulation at this node is perhaps over-parameterized?\n",
                                                        "## Fitting the glmm at this node using glmer() in lme4 may provide more information\n",
                                                        "## If glmer() can fit this model then please send a bug report to fraseriain.lewis@uzh.ch\n",sep="");
                                                         stop("");}

                                   used.inla<-FALSE;## flip
                                    
                                 } else {if (verbose) cat("Using INLA (glmm)\n");
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
                            
                            if(used.inla==FALSE){## organize output from C
                            res.list[[child.name]]<-res.c[[1]][1];            
                            mymodes[[child]]<-res.c[[1]][-c(1:3)];## remove mlik - this is first entry, and error code and hessian accuracy
                            mymodes[[child]]<-mymodes[[child]][which(mymodes[[child]]!=.Machine$double.xmax)];## this discards all "empty" parameters 
                            error.code[child]<-res.c[[1]][2];
                            hessian.accuracy[child]<-res.c[[1]][3];
                            INLA.marginals<-c(INLA.marginals,FALSE);
                            } else {
                                    ## organize output from INLA
                                    res.list[[child.name]]<-res.inla$mlik[2];## [2] is for Gaussian rather than Integrated estimate 
                                    mymodes[[child]]<-getModeVector(list.fixed=res.inla$marginals.fixed,list.hyper=res.inla$marginals.hyperpar);
                                    mymargs[[child]]<-getMargsINLA(list.fixed=res.inla$marginals.fixed,list.hyper=res.inla$marginals.hyperpar);
                                    #print(mymodes);stop("");
                                    error.code[child]<-NA;## not available from INLA
                                    hessian.accuracy[child]<-NA;## not available from INLA
                                    INLA.marginals<-c(INLA.marginals,TRUE);
                            } 
                            
                            nom<-colnames(dag.m)[which(dag.m[child,]==1)];
                            if(var.types[child]=="1"){nom<-c("(Intercept)",nom,"group.precision");}## binomial : just some naming for use later
                            if(var.types[child]=="2" && !(child%in%grouped.vars)){nom<-c("(Intercept)",nom,"precision","precision");} ## gaus and not grouped
                            if(var.types[child]=="2" && (child%in%grouped.vars)){nom<-c("(Intercept)",nom,"group.precision","precision");} 
                            if(var.types[child]=="3"){nom<-c("(Intercept)",nom,"group.precision");} ## pois}
                            mynom<-NULL;
                            for(j in 1:length(mymodes[[child]])){mynom<-c(mynom,paste(colnames(dag.m)[child],nom[j],sep="|"));
                                                                 }
                            names(mymodes[[child]])<-mynom; 
                            if(used.inla==TRUE){names(mymargs[[child]])<-mynom;}
                            
                            if(!is.null(orig.force.method)){force.method<-orig.force.method;} ## reset force.method after INLA crash
                            ############################################################
                            ## Finished with current node
                            ############################################################

 } ## end of nodes loop
                           
      #################################################
      ## Further tidy up of results across all nodes
      #################################################
      names(mymodes)<-colnames(dag.m);
      
      res.list[["modes"]]<-mymodes;
      res.list[["error.code"]]<-error.code;
      res.list[["hessian.accuracy"]]<-hessian.accuracy;
      names(res.list[["error.code"]])<-names(res.list[["hessian.accuracy"]])<-colnames(dag.m);
      res.list[["error.code.desc"]]<-as.character(res.list[["error.code"]]);
      res.list[["error.code.desc"]]<-ifelse(res.list[["error.code.desc"]]==0,"success",res.list[["error.code"]]);
      res.list[["error.code.desc"]]<-ifelse(res.list[["error.code.desc"]]==1,"warning: mode results may be unreliable (optimiser terminated unusually)",res.list[["error.code.desc"]]);
      res.list[["error.code.desc"]]<-ifelse(res.list[["error.code.desc"]]==2,"error - logscore is NA - model could not be fitted",res.list[["error.code.desc"]]);
      res.list[["error.code.desc"]]<-ifelse(res.list[["error.code.desc"]]==4,"warning: fin.diff hessian estimation terminated unusually ",res.list[["error.code.desc"]]);
      res.list[["mlik"]]<-sum(unlist(res.list[1:dim(dag.m)[1]]));## overall mlik 

      res.list[["used.INLA"]]<-INLA.marginals;## vector - TRUE if INLA used false otherwise

      #####
      ##### Additional part only run if user wants marginal distributions
      #####
      if(compute.fixed){ if (verbose) cat("Processing marginal distributions for non-INLA nodes...\n");
                         ## might have some already computed from INLA
                         if(length(which(INLA.marginals==FALSE))==0){## ALL INLA so finished
                            names(mymargs)<-colnames(dag.m)[which(INLA.marginals==TRUE)];
                            res.list[["marginals"]]<-mymargs;
                         } else {## At least one C and this creates its own res.list[["marginals"]]
                         ## now get rest using C
                         max.parents<-max(apply(dag.m,1,sum));## over all nodes - different from above
                         res.list<-getmarginals(res.list, ## rest of arguments as for call to C fitabn
                                                data.df,
                                                dag.m,
                                                var.types,
                                                max.parents,
                                                mean,prec,loggam.shape,loggam.inv.scale,
                                                max.iters,epsabs,verbose,error.verbose,
                                                as.integer(grouped.vars-1),## int.vector of variables which are mixed model nodes -1 for C (changed from earlier fitabn)
                                                group.ids,
                                                epsabs.inner,max.iters.inner,finite.step.size,hessian.params,max.iters.hessian,
                                                min.pdf,marginal.node,marginal.param,variate.vec,n.grid,
                                                INLA.marginals,
                                                max.grid.iter,
                                               as.double(max.hessian.error),
                                               as.double(factor.brent),
                                               as.integer(maxiters.hessian.brent),
                                               as.double(num.intervals.brent) 
                                                );
                                     }
                        
                        ## at least one INLA node so we need to combine 
                        ## res.list[["inla.margs"]] with res.list[["marginals"]] 
                        if(length(which(INLA.marginals==TRUE))>0){
                            names(mymargs)<-colnames(dag.m)[which(INLA.marginals==TRUE)];
                            res.list[["inla.margs"]]<-mymargs;
                            ## also have res.list[["marginals"]] from getmarginalsC above 
                           
                            masterlist<-list();
                            for(i in rownames(dag.m)){
                                     attempt1<-which(names(res.list[["inla.margs"]])==i);
                                     attempt2<-which(names(res.list[["marginals"]])==i);
                                     if(length(attempt1)>0){
                                                            masterlist[[i]]<-res.list[["inla.margs"]][[attempt1]];
                                     } else {masterlist[[i]]<-res.list[["marginals"]][[attempt2]];}
                            #print(masterlist);
                            }
                            
                            #res.list[["marginals.master"]]<-masterlist;
                            res.list[["marginals"]]<-masterlist;
                            res.list[["inla.margs"]]<-NULL;

                        }

        ## Three final optional operations
        ## 1. evaluate density across an equally spaced grid - this used spline interpolation
        ## 2. standardise the area to unity - this should alread be very close but will do no harm
        ## 3. compute quantiles - this is done after 1. and 2. (assuming they were turned on
        
        ## 1.               
        if(!is.null(n.grid)){## evaluated density across an equal grid - use spline interpolation
                             res.list[["marginals"]]<-eval.across.grid(res.list[["marginals"]],n.grid,marginal.node);
                             }
        ## 2.
        if(std.area){## want to standardize area under curve to unity - might be slightly adrift as is depending on accuracy of approx's
                    if(is.null(n.grid)){stop("must provide n.grid if using std.area!");}
                    res.list[["marginals"]]<-std.area.under.grid(res.list[["marginals"]],marginal.node);
                     } 
        ## 3.
        if(!is.null(marginal.quantiles)){res.list[["marginal.quantiles"]]<-get.quantiles(res.list[["marginals"]],marginal.quantiles,marginal.node);}
 
        } ## end of compute.fixed 

      #########################################################
      ## Rgraph/graphviz part
      #########################################################
      if(create.graph){
      if(!requireNamespace("Rgraphviz", quietly=TRUE)){stop("library Rgraphviz is not available!\nRgraphviz is available as part of the bioconductor project - see http://www.bioconductor.org/install\nRgraphviz is required is create.graph=TRUE");}
      
      mygraph<-new("graphAM",adjMat=t(dag.m),edgemode="directed");
      res.list[["graph"]]<-mygraph;

      }


if (verbose) cat("########End of DAG fitting #############################\n");
if (verbose) cat("########################################################\n");
return(res.list);

}

#####################################################
## function to extract the mode from INLA output
#####################################################
getModeVector<-function(list.fixed,list.hyper){
  
 ## listfixed is a list of matrices of two cols x, y
 modes<-NULL;
 for(i in 1:length(list.fixed)){
   mymarg<-list.fixed[[i]];## matrix 2 cols
   mymarg<-spline(mymarg);
   modes<-c(modes,mymarg$x[which(mymarg$y==max(mymarg$y))]);## find x value which corresponds to the max y value
 }
 
 ## now for hyperparam if there is one
 if(is.list(list.hyper)){## might be no hyperparam e.g. binomial or poisson glmm
      if(length(list.hyper)==1){
      mymarg<-list.hyper[[1]];## matrix 2 cols
      mymarg<-spline(mymarg);
      modes<-c(modes,mymarg$x[which(mymarg$y==max(mymarg$y))]);}## find x value which corresponds to the max y value
     
       if(length(list.hyper)==2){ ## gaussian glm
      mymarg<-list.hyper[[2]];## matrix 2 cols - group level precision
      mymarg<-spline(mymarg);
      modes<-c(modes,mymarg$x[which(mymarg$y==max(mymarg$y))]); 
      mymarg<-list.hyper[[1]];## matrix 2 cols - residual level precision
      mymarg<-spline(mymarg);
      modes<-c(modes,mymarg$x[which(mymarg$y==max(mymarg$y))]);
      ### [[2]] then [[1]] to keep same order a C code
       }
}

return(modes);

}

#####################################################
## function to extract the mode from INLA output
#####################################################
getMargsINLA<-function(list.fixed,list.hyper){
  
 ## listfixed is a list of matrices of two cols x, y
 margs<-list();
 for(i in 1:length(list.fixed)){
   margs[[i]]<-list.fixed[[i]];## matrix 2 cols
 }  
 ## now for hyperparam if there is one
 if(is.list(list.hyper)){## might be no hyperparam e.g. binomial or poisson glmm
    if(length(list.hyper)==1){## poisson or binomial glmm 
    margs[[i+1]]<-list.hyper[[1]];}## matrix 2 cols
    if(length(list.hyper)==2){## gaussian glmm two entries
    margs[[i+1]]<-list.hyper[[2]]; ## this is the group level precision 
    margs[[i+2]]<-list.hyper[[1]];}## this is the residual precision
                   ## note this is [[2]] then [[1]] to keep same order a C code
     
 }

return(margs);

}
#####################################################
## function to get marginal across an equal grid
#####################################################
eval.across.grid<-function(mylist,n.grid,single){

if(is.null(single)){
for(i in 1:length(mylist)){## for each node
       q.inner.list<-mylist[[i]];##copy
       for(j in 1:length(q.inner.list)){## for each parameter
             mymat<-q.inner.list[[j]];## copy - this is a matrix
             q.mat<-matrix(data=rep(NA,2*n.grid),ncol=2);## create new matrix
             colnames(q.mat)<-c("x","f(x)");
             interp<-spline(mymat,n=n.grid);## interpolate over equal spaced grid of n points
             q.mat[,1]<-interp$x;
             q.mat[,2]<-interp$y;
             q.inner.list[[j]]<-q.mat;## overwrite 
         }
         mylist[[i]]<-q.inner.list;
}
} else {## have only a single node and parameter
             mymat<-mylist[[1]];
             q.mat<-matrix(data=rep(NA,2*n.grid),ncol=2);## create new matrix
             colnames(q.mat)<-c("x","f(x)");
             interp<-spline(mymat,n=n.grid);## interpolate over equal spaced grid of n points
             q.mat[,1]<-interp$x;
             q.mat[,2]<-interp$y;
             mylist[[1]]<-q.mat;## overwrite 
} 



return(mylist);
}
##################################################################################
## function to get std. are under marginal to exactly unity
## it should be very close to unity but in some cases due to numerical accuracy
## differences (since each point is a separate estimate) this might be a little adrift
## turn this option off to see how reliable the original estimation is
##################################################################################
std.area.under.grid<-function(mylist,single){
 
if(is.null(single)){        
for(i in 1:length(mylist)){## for each node
       q.inner.list<-mylist[[i]];##copy
       for(j in 1:length(q.inner.list)){## for each parameter
             mymat<-q.inner.list[[j]];## copy - this is a matrix
             cur.area<-(mymat[2,1]-mymat[1,1])*sum(mymat[,2]);## delta.x used - equal sized grid
             mymat[,2]<-mymat[,2]/cur.area;## std to ~ 1.0
             q.inner.list[[j]]<-mymat;## overwrite 
         }
         mylist[[i]]<-q.inner.list;
}
} else {
             mymat<-mylist[[1]];## copy - this is a matrix
             cur.area<-(mymat[2,1]-mymat[1,1])*sum(mymat[,2]);## delta.x used - equal sized grid
             mymat[,2]<-mymat[,2]/cur.area;## std to ~ 1.0
             mylist[[1]]<-mymat;## overwrite 

         }

return(mylist);
}
#####################################################
## function to quantiles
#####################################################
get.quantiles<-function(mylist,quantiles, single){

if(is.null(single)){  
 for(i in 1:length(mylist)){
       q.inner.list<-mylist[[i]];##copy
       for(j in 1:length(q.inner.list)){
             mymat<-q.inner.list[[j]];## copy - this is a matrix
             q.mat<-matrix(data=rep(NA,2*length(quantiles)),ncol=2);## create new matrix
             colnames(q.mat)<-c("P(X<=x)","x");
             q.mat[,1]<-quantiles;
             q.mat<-get.ind.quantiles(q.mat,mymat);## the actual quantile arithmetic
             q.inner.list[[j]]<-q.mat;## overwrite 
         }
         mylist[[i]]<-q.inner.list;
}
} else { 
             mymat<-mylist[[1]];## copy - this is a matrix
             q.mat<-matrix(data=rep(NA,2*length(quantiles)),ncol=2);## create new matrix
             colnames(q.mat)<-c("P(X<=x)","x");
             q.mat[,1]<-quantiles;
             q.mat<-get.ind.quantiles(q.mat,mymat);## the actual quantile arithmetic
             mylist[[1]]<-q.mat;## overwrite 

  }
  return(mylist);           
}

## helper function for get.quantiles above
get.ind.quantiles<-function(outmat,inmat){
      ##outmat is a matrix where the first col has the desired quantiles
      ## we want to estimate this and out in into the second col
      ## inmat is the actual x,f(x) matrix
      qs.vec<-outmat[,1];
      x<-inmat[,1];
      fx<-inmat[,2];
      cum.den<-cumsum(fx);## cumulative of f(x) values
      sum.den<-sum(fx);## total sum of f(x)
      cum.f<-cum.den/sum.den;## cumulative density function
      row<-1;
      for(qs in qs.vec){## for each quantile
       ## find row in inmat which has cumulative f(x)>q.s
       outmat[row,2]<-x[which(cum.f>qs)[1]];## find first row to exceed quantile value q.s and get x value
       row<-row+1;
      }

  return(outmat);
}
