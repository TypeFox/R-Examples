#################################################################################################
## getmarginals.R
## function for computing marginal posterior densities using C and is called from fit.dag()
## Only to be called internally.
#################################################################################################

getmarginals<-function(res.list, ## rest of arguments as for call to C fitabn
                       data.df,dag.m,var.types,max.parents,
                       mean,prec,loggam.shape,loggam.inv.scale,
                       max.iters,epsabs,verbose,error.verbose,
                       grouped.vars,## int.vector of variables which are mixed model nodes
                       group.ids,
                       epsabs.inner,max.iters.inner,finite.step.size,hessian.params,max.iters.hessian,min.pdf,marginal.node,marginal.param, variate.vec,n.grid,
                       INLA.marginals, ## this argument contains the marginals already computed via calls to INLA
                       iter.max, ## max. number of steps to take in one direction during evaluation - to avoid taking forever
                       max.hessian.error,
                       factor.brent,
                       maxiters.hessian.brent,
                       num.intervals.brent){

if(!is.null(marginal.node) && is.null(variate.vec)){stop("must supply variate.vec if using a single node!");}

marginals<-list();

if( !is.null(marginal.node)) {## in single node case 
    if(res.list[["error.code"]][marginal.node]!=0){stop("---- Cannot compute marginal density as the mlik is unreliable for this node.\nTry re-running with different finite difference parameters or else increase max.hessian.error ----");}
    nodeid<-marginal.node;
    paramid<-marginal.param;
    if (verbose) cat("processing ",names(res.list$modes[[nodeid]])[paramid],"\n");
    curnom<-names(res.list$modes[[nodeid]])[paramid];
    first<-TRUE;
    for(betafixed in variate.vec){
               marg.res <- .Call("fitabn_marginals",data.df,as.integer(dag.m),as.integer(dim(data.df)[2]),
                            as.integer(var.types),as.integer(max.parents),
                            as.double(mean),as.double(1/sqrt(prec)),as.double(loggam.shape),as.double(1/loggam.inv.scale),
                            as.integer(max.iters),as.double(epsabs),
                            as.integer(verbose),as.integer(error.verbose),
                            grouped.vars,## int.vector of variables which are mixed model nodes
                            group.ids,
                            as.double(epsabs.inner),
                            as.integer(max.iters.inner),
                            as.double(finite.step.size),
                            as.double(hessian.params),
                            ## additional marginal arguments - node, parameter, and relevant modes
                            as.integer(nodeid-1),## C index from 0
                            as.integer(paramid-1),## C index from 0
                            as.double(res.list$modes[[nodeid]]), ## relevant mode estimates
                            as.double(betafixed), ## value to be evaluated e.g. want f(x) given x - betafixed is x
                            as.double(res.list[[names(res.list$modes)[nodeid]]]),## mlik for the node
                            as.integer(max.iters.hessian), 
                            as.double(max.hessian.error),
                            as.double(factor.brent),
                            as.integer(maxiters.hessian.brent),
                            as.double(num.intervals.brent),

              PACKAGE="abn" ## uncomment to load as package not shlib
              );
          if(first){tmp<-matrix(data=c(betafixed,marg.res),ncol=2,byrow=TRUE);colnames(tmp)<-c("x","f(x)");
                    marginals[[curnom]]<-tmp;first<-FALSE;
          } else{ marginals[[curnom]]<-rbind(marginals[[curnom]],c(betafixed,marg.res));}
          if (verbose) cat("evaluating at x=",betafixed," f(x)=",marg.res,"\n",sep="");
          }
         marginals[[curnom]]<-marginals[[curnom]][order(marginals[[curnom]][,1]),];

} else {

  ## for each node
  index<-1;
  for(nodeid in 1:length(res.list$modes)){
   
  # nodeid=6;
    if(!INLA.marginals[nodeid]){## COMPUTE NODE MARGINALS if NOT already done via INLA
     node.marginals<-list();
     
     ## for each parameter
     for(paramid in 1:length(res.list$modes[[nodeid]])){

   #paramid<-3;    
         if (verbose) cat("processing ",names(res.list$modes[[nodeid]])[paramid],"\n");      
         curnom<-names(res.list$modes[[nodeid]])[paramid]; 
         have.precision<-ifelse(grepl("prec",curnom),TRUE,FALSE);
         if(res.list[["error.code"]][nodeid]!=0){ if (verbose) cat("--- NOTE DROPPING parameter: ",curnom," as mlik is unreliable for this node.\nTry re-running with different finite difference parameters---\n");}  
         betafixedMode<-res.list$modes[[nodeid]][paramid];## just evaluate at the mode for paramid in node id
         #betafixedMode<-signif(betafixedMode,1);## adjust the mode a little as this seems to cause problems in the C optim in some unusual cases
         ##############################################################################################################
         ## STEP 1. get the value at the mode and one point either side (n.b. this may not give a turning point....)
         ## GET INITIAL POINTS EITHER SIDE OF MODE
         ## README THIS NEXT SECTION IS LONG BUT JUST REPETITION!
         ##############################################################################################################
         ## this is crucial as it controls the rest - we start with a guess and then if this is too wide or two small then change
         if(abs(betafixedMode)<1E-05){## if have a centred variable so create two evaluation point either side using se.delta 
                                      x.delta<-se.delta<-(sd(data.df[,nodeid])/sqrt(dim(data.df)[1]))/10;## a feeble effort to get a good step size
                                      variate.x<-c(betafixedMode-x.delta,betafixedMode,betafixedMode+x.delta);
         } else {## not centred variable to again try one point either side of mode
                                      x.delta<-betafixedMode*0.025
                                      variate.x<-c(betafixedMode-x.delta,betafixedMode,betafixedMode+x.delta);
                                      if(have.precision && variate.x[1]<0.001){variate.x[1]<-0.002;}## do not evaluate precision at negative x!
          }
          ###############################################################################################################
          ### This is the first step is a long section of repetition which are just successive attempts to find a good min. step size
          ### for moving across the x grid. Currently 10 attempts each successive *2 or /2 - final one used if nothing else works
          ###############################################################################################################
          ## FIRST TRY
          ## created initial guess x points now se how goo they are and adjust if necessary
          row<-1;
          mymat<-matrix(data=c(variate.x,rep(NA,3)),ncol=2,byrow=TRUE);colnames(mymat)<-c("x","f(x)");## storage for three rows
          for(betafixed in variate.x){## for each of the three initial points mode +/- one other
          
          marg.res <- .Call("fitabn_marginals",data.df,as.integer(dag.m),as.integer(dim(data.df)[2]),
                            as.integer(var.types),as.integer(max.parents),
                            as.double(mean),as.double(1/sqrt(prec)),as.double(loggam.shape),as.double(1/loggam.inv.scale),
                            as.integer(max.iters),as.double(epsabs),
                            as.integer(verbose),as.integer(error.verbose),
                            grouped.vars,## int.vector of variables which are mixed model nodes
                            group.ids,
                            as.double(epsabs.inner),
                            as.integer(max.iters.inner),
                            as.double(finite.step.size),
                            as.double(hessian.params),
                            ## additional marginal arguments - node, parameter, and relevant modes
                            as.integer(nodeid-1),## C index from 0
                            as.integer(paramid-1),## C index from 0
                            as.double(res.list$modes[[nodeid]]), ## relevant mode estimates
                            as.double(betafixed), ## value to be evaluated e.g. want f(x) given x - betafixed is x
                            as.double(res.list[[names(res.list$modes)[nodeid]]]),## mlik for the node
                            as.integer(max.iters.hessian),as.double(max.hessian.error),
                            as.double(factor.brent),
                            as.integer(maxiters.hessian.brent),
                            as.double(num.intervals.brent),
              PACKAGE="abn" ## uncomment to load as package not shlib
              );
              #cat("betafixed=",betafixed," gvalue=",marg.res,"\n");
              mymat[row,]<-c(betafixed,marg.res);row<-row+1; 
              }
              #cat("original points\n");print(mymat);
              gvalue<-mymat[,2];## f(x) values
              big.ok<-small.ok<-FALSE;
              #cat("diffs:", abs(gvalue[1]-gvalue[2])," ",abs(gvalue[3]-gvalue[2])," 0.02*gvalue=",0.02*gvalue[2]," 0.1*gvalue=",0.1*gvalue[2],"\n");
              #stop("");
              if(   abs(gvalue[1]-gvalue[2])> 0.02*gvalue[2] && abs(gvalue[3]-gvalue[2])> 0.02*gvalue[2] ){big.ok<-TRUE;}
              if(   abs(gvalue[1]-gvalue[2])< 0.1*gvalue[2] && abs(gvalue[3]-gvalue[2])< 0.1*gvalue[2] ){small.ok<-TRUE;}## ok choices of variate.x
              if(!big.ok){#cat("increasing x.delta by *2 from ",x.delta," ");
                          x.delta<-x.delta*2;#cat("to ",x.delta,"\n");
                         }## gap is not big enough so increase x.delta by factor of 5
              if(!small.ok){
                           #cat("decreasing x.delta by *2 from ",x.delta," ");
                          x.delta<-x.delta/2;#cat("to ",x.delta,"\n");
                            } ## gap too big so decrease x.delta by factor of 5              
              
              ## SECOND try for x.delta
              if( !(big.ok && small.ok)){
              variate.x<-c(betafixedMode-x.delta,betafixedMode,betafixedMode+x.delta);
              row<-1;
              mymat<-matrix(data=c(variate.x,rep(NA,3)),ncol=2,byrow=TRUE);colnames(mymat)<-c("x","f(x)");## storage for three rows
              for(betafixed in variate.x){## for each of the three initial points mode +/- one other
              gvalue <- .Call("fitabn_marginals",data.df,as.integer(dag.m),as.integer(dim(data.df)[2]),
                            as.integer(var.types),as.integer(max.parents),
                            as.double(mean),as.double(1/sqrt(prec)),as.double(loggam.shape),as.double(1/loggam.inv.scale),
                            as.integer(max.iters),as.double(epsabs),
                            as.integer(verbose),as.integer(error.verbose),
                            grouped.vars,## int.vector of variables which are mixed model nodes
                            group.ids,
                            as.double(epsabs.inner),
                            as.integer(max.iters.inner),
                            as.double(finite.step.size),
                            as.double(hessian.params),
                            ## additional marginal arguments - node, parameter, and relevant modes
                            as.integer(nodeid-1),## C index from 0
                            as.integer(paramid-1),## C index from 0
                            as.double(res.list$modes[[nodeid]]), ## relevant mode estimates
                            as.double(betafixed), ## value to be evaluated e.g. want f(x) given x - betafixed is x
                            as.double(res.list[[names(res.list$modes)[nodeid]]]),## mlik for the node
                            as.integer(max.iters.hessian),as.double(max.hessian.error),
                            as.double(factor.brent),
                            as.integer(maxiters.hessian.brent),
                            as.double(num.intervals.brent),
                  PACKAGE="abn" ## uncomment to load as package not shlib
                  );
                  mymat[row,]<-c(betafixed,gvalue);row<-row+1; 
                  }
                  gvalue<-mymat[,2];## f(x) values
                  big.ok<-small.ok<-FALSE;
                  #cat("diffs:", abs(gvalue[1]-gvalue[2])," ",abs(gvalue[3]-gvalue[2])," 0.02*gvalue=",0.02*gvalue[2]," 0.1*gvalue=",0.1*gvalue[2],"\n");
                  if(   abs(gvalue[1]-gvalue[2])> 0.02*gvalue[2] && abs(gvalue[3]-gvalue[2])> 0.02*gvalue[2] ){big.ok<-TRUE;}
                  if(   abs(gvalue[1]-gvalue[2])< 0.1*gvalue[2] && abs(gvalue[3]-gvalue[2])< 0.1*gvalue[2] ){small.ok<-TRUE;}## ok choices of variate.x
                  if(!big.ok){
                              #cat("increasing again x.delta by *2 from ",x.delta," ");
                              x.delta<-x.delta*2;#cat("to ",x.delta,"\n");
                              }## gap is not big enough so increase x.delta by factor of 5
                   if(!small.ok){
                              #cat("decreasing again x.delta by *2 from ",x.delta," ");
                              x.delta<-x.delta/2;#cat("to ",x.delta,"\n");
                               } ## gap too big so decrease x.delta by factor of 5              

              }
              
              ## THIRD try for x.delta
              if( !(big.ok && small.ok)){
              variate.x<-c(betafixedMode-x.delta,betafixedMode,betafixedMode+x.delta);
              row<-1;
              mymat<-matrix(data=c(variate.x,rep(NA,3)),ncol=2,byrow=TRUE);colnames(mymat)<-c("x","f(x)");## storage for three rows
              for(betafixed in variate.x){## for each of the three initial points mode +/- one other
              gvalue <- .Call("fitabn_marginals",data.df,as.integer(dag.m),as.integer(dim(data.df)[2]),
                            as.integer(var.types),as.integer(max.parents),
                            as.double(mean),as.double(1/sqrt(prec)),as.double(loggam.shape),as.double(1/loggam.inv.scale),
                            as.integer(max.iters),as.double(epsabs),
                            as.integer(verbose),as.integer(error.verbose),
                            grouped.vars,## int.vector of variables which are mixed model nodes
                            group.ids,
                            as.double(epsabs.inner),
                            as.integer(max.iters.inner),
                            as.double(finite.step.size),
                            as.double(hessian.params),
                            ## additional marginal arguments - node, parameter, and relevant modes
                            as.integer(nodeid-1),## C index from 0
                            as.integer(paramid-1),## C index from 0
                            as.double(res.list$modes[[nodeid]]), ## relevant mode estimates
                            as.double(betafixed), ## value to be evaluated e.g. want f(x) given x - betafixed is x
                            as.double(res.list[[names(res.list$modes)[nodeid]]]),## mlik for the node
                            as.integer(max.iters.hessian),as.double(max.hessian.error),
                            as.double(factor.brent),
                            as.integer(maxiters.hessian.brent),
                            as.double(num.intervals.brent),
                  PACKAGE="abn" ## uncomment to load as package not shlib
                  );
                  mymat[row,]<-c(betafixed,gvalue);row<-row+1; 
                  }
                  gvalue<-mymat[,2];## f(x) values
                  big.ok<-small.ok<-FALSE;
                  #cat("diffs:", abs(gvalue[1]-gvalue[2])," ",abs(gvalue[3]-gvalue[2])," 0.02*gvalue=",0.02*gvalue[2]," 0.1*gvalue=",0.1*gvalue[2],"\n");
                  if(   abs(gvalue[1]-gvalue[2])> 0.02*gvalue[2] && abs(gvalue[3]-gvalue[2])> 0.02*gvalue[2] ){big.ok<-TRUE;}
                  if(   abs(gvalue[1]-gvalue[2])< 0.1*gvalue[2] && abs(gvalue[3]-gvalue[2])< 0.1*gvalue[2] ){small.ok<-TRUE;}## ok choices of variate.x
                  if(!big.ok){
                              #cat("increasing 3rd time x.delta by *2 from ",x.delta," ");
                              x.delta<-x.delta*2;#cat("to ",x.delta,"\n");
                             }## gap is not big enough so increase x.delta by factor of 5
                   if(!small.ok){
                              #cat("decreasing 3rd time x.delta by *2 from ",x.delta," ");
                              x.delta<-x.delta/2;#cat("to ",x.delta,"\n");
                                } ## gap too big so decrease x.delta by factor of 5              

              }
 
              ## FOURTH try for x.delta
              if( !(big.ok && small.ok)){
              variate.x<-c(betafixedMode-x.delta,betafixedMode,betafixedMode+x.delta);
              row<-1;
              mymat<-matrix(data=c(variate.x,rep(NA,3)),ncol=2,byrow=TRUE);colnames(mymat)<-c("x","f(x)");## storage for three rows
              for(betafixed in variate.x){## for each of the three initial points mode +/- one other
              gvalue <- .Call("fitabn_marginals",data.df,as.integer(dag.m),as.integer(dim(data.df)[2]),
                            as.integer(var.types),as.integer(max.parents),
                            as.double(mean),as.double(1/sqrt(prec)),as.double(loggam.shape),as.double(1/loggam.inv.scale),
                            as.integer(max.iters),as.double(epsabs),
                            as.integer(verbose),as.integer(error.verbose),
                            grouped.vars,## int.vector of variables which are mixed model nodes
                            group.ids,
                            as.double(epsabs.inner),
                            as.integer(max.iters.inner),
                            as.double(finite.step.size),
                            as.double(hessian.params),
                            ## additional marginal arguments - node, parameter, and relevant modes
                            as.integer(nodeid-1),## C index from 0
                            as.integer(paramid-1),## C index from 0
                            as.double(res.list$modes[[nodeid]]), ## relevant mode estimates
                            as.double(betafixed), ## value to be evaluated e.g. want f(x) given x - betafixed is x
                            as.double(res.list[[names(res.list$modes)[nodeid]]]),## mlik for the node
                            as.integer(max.iters.hessian),as.double(max.hessian.error),
                            as.double(factor.brent),
                            as.integer(maxiters.hessian.brent),
                            as.double(num.intervals.brent),
                  PACKAGE="abn" ## uncomment to load as package not shlib
                  );
                  mymat[row,]<-c(betafixed,gvalue);row<-row+1; 
                  }
                  gvalue<-mymat[,2];## f(x) values
                  big.ok<-small.ok<-FALSE;
                  #cat("diffs:", abs(gvalue[1]-gvalue[2])," ",abs(gvalue[3]-gvalue[2])," 0.02*gvalue=",0.02*gvalue[2]," 0.1*gvalue=",0.1*gvalue[2],"\n");
                  if(   abs(gvalue[1]-gvalue[2])> 0.02*gvalue[2] && abs(gvalue[3]-gvalue[2])> 0.02*gvalue[2] ){big.ok<-TRUE;}
                  if(   abs(gvalue[1]-gvalue[2])< 0.1*gvalue[2] && abs(gvalue[3]-gvalue[2])< 0.1*gvalue[2] ){small.ok<-TRUE;}## ok choices of variate.x
                  if(!big.ok){
                              #cat("increasing 4th time x.delta by *2 from ",x.delta," ");
                              x.delta<-x.delta*2;#cat("to ",x.delta,"\n");
                              }## gap is not big enough so increase x.delta by factor of 5
                   if(!small.ok){
                              #cat("decreasing 4th time x.delta by *2 from ",x.delta," ");
                              x.delta<-x.delta/2;#cat("to ",x.delta,"\n");
                                } ## gap too big so decrease x.delta by factor of 5              

              }

              ## FIFTH try for x.delta
              if( !(big.ok && small.ok)){
              variate.x<-c(betafixedMode-x.delta,betafixedMode,betafixedMode+x.delta);
              row<-1;
              mymat<-matrix(data=c(variate.x,rep(NA,3)),ncol=2,byrow=TRUE);colnames(mymat)<-c("x","f(x)");## storage for three rows
              for(betafixed in variate.x){## for each of the three initial points mode +/- one other
              gvalue <- .Call("fitabn_marginals",data.df,as.integer(dag.m),as.integer(dim(data.df)[2]),
                            as.integer(var.types),as.integer(max.parents),
                            as.double(mean),as.double(1/sqrt(prec)),as.double(loggam.shape),as.double(1/loggam.inv.scale),
                            as.integer(max.iters),as.double(epsabs),
                            as.integer(verbose),as.integer(error.verbose),
                            grouped.vars,## int.vector of variables which are mixed model nodes
                            group.ids,
                            as.double(epsabs.inner),
                            as.integer(max.iters.inner),
                            as.double(finite.step.size),
                            as.double(hessian.params),
                            ## additional marginal arguments - node, parameter, and relevant modes
                            as.integer(nodeid-1),## C index from 0
                            as.integer(paramid-1),## C index from 0
                            as.double(res.list$modes[[nodeid]]), ## relevant mode estimates
                            as.double(betafixed), ## value to be evaluated e.g. want f(x) given x - betafixed is x
                            as.double(res.list[[names(res.list$modes)[nodeid]]]),## mlik for the node
                            as.integer(max.iters.hessian),as.double(max.hessian.error),
                            as.double(factor.brent),
                            as.integer(maxiters.hessian.brent),
                            as.double(num.intervals.brent),
                  PACKAGE="abn" ## uncomment to load as package not shlib
                  );
                  mymat[row,]<-c(betafixed,gvalue);row<-row+1; 
                  }
                  gvalue<-mymat[,2];## f(x) values
                  big.ok<-small.ok<-FALSE;
                  #cat("diffs:", abs(gvalue[1]-gvalue[2])," ",abs(gvalue[3]-gvalue[2])," 0.02*gvalue=",0.02*gvalue[2]," 0.1*gvalue=",0.1*gvalue[2],"\n");
                  if(   abs(gvalue[1]-gvalue[2])> 0.02*gvalue[2] && abs(gvalue[3]-gvalue[2])> 0.02*gvalue[2] ){big.ok<-TRUE;}
                  if(   abs(gvalue[1]-gvalue[2])< 0.1*gvalue[2] && abs(gvalue[3]-gvalue[2])< 0.1*gvalue[2] ){small.ok<-TRUE;}## ok choices of variate.x
                  if(!big.ok){
                              #cat("increasing 5th time x.delta by *2 from ",x.delta," ");
                              x.delta<-x.delta*2;#cat("to ",x.delta,"\n");
                             }## gap is not big enough so increase x.delta by factor of 5
                   if(!small.ok){
                              #cat("decreasing 5th time x.delta by *2 from ",x.delta," ");
                              x.delta<-x.delta/2;#cat("to ",x.delta,"\n");
                                 } ## gap too big so decrease x.delta by factor of 5              

              }
              ## SIXTH try for x.delta
              if( !(big.ok && small.ok)){
              variate.x<-c(betafixedMode-x.delta,betafixedMode,betafixedMode+x.delta);
              row<-1;
              mymat<-matrix(data=c(variate.x,rep(NA,3)),ncol=2,byrow=TRUE);colnames(mymat)<-c("x","f(x)");## storage for three rows
              for(betafixed in variate.x){## for each of the three initial points mode +/- one other
              gvalue <- .Call("fitabn_marginals",data.df,as.integer(dag.m),as.integer(dim(data.df)[2]),
                            as.integer(var.types),as.integer(max.parents),
                            as.double(mean),as.double(1/sqrt(prec)),as.double(loggam.shape),as.double(1/loggam.inv.scale),
                            as.integer(max.iters),as.double(epsabs),
                            as.integer(verbose),as.integer(error.verbose),
                            grouped.vars,## int.vector of variables which are mixed model nodes
                            group.ids,
                            as.double(epsabs.inner),
                            as.integer(max.iters.inner),
                            as.double(finite.step.size),
                            as.double(hessian.params),
                            ## additional marginal arguments - node, parameter, and relevant modes
                            as.integer(nodeid-1),## C index from 0
                            as.integer(paramid-1),## C index from 0
                            as.double(res.list$modes[[nodeid]]), ## relevant mode estimates
                            as.double(betafixed), ## value to be evaluated e.g. want f(x) given x - betafixed is x
                            as.double(res.list[[names(res.list$modes)[nodeid]]]),## mlik for the node
                            as.integer(max.iters.hessian),as.double(max.hessian.error),
                            as.double(factor.brent),
                            as.integer(maxiters.hessian.brent),
                            as.double(num.intervals.brent),
                  PACKAGE="abn" ## uncomment to load as package not shlib
                  );
                  mymat[row,]<-c(betafixed,gvalue);row<-row+1; 
                  }
                  gvalue<-mymat[,2];## f(x) values
                  big.ok<-small.ok<-FALSE;
                  #cat("diffs:", abs(gvalue[1]-gvalue[2])," ",abs(gvalue[3]-gvalue[2])," 0.02*gvalue=",0.02*gvalue[2]," 0.1*gvalue=",0.1*gvalue[2],"\n");
                  if(   abs(gvalue[1]-gvalue[2])> 0.02*gvalue[2] && abs(gvalue[3]-gvalue[2])> 0.02*gvalue[2] ){big.ok<-TRUE;}
                  if(   abs(gvalue[1]-gvalue[2])< 0.1*gvalue[2] && abs(gvalue[3]-gvalue[2])< 0.1*gvalue[2] ){small.ok<-TRUE;}## ok choices of variate.x
                  if(!big.ok){
                              #cat("increasing 6th time x.delta by *2 from ",x.delta," ");
                              x.delta<-x.delta*2;#cat("to ",x.delta,"\n");
                             }## gap is not big enough so increase x.delta by factor of 5
                   if(!small.ok){
                              #cat("decreasing 6th time x.delta by *2 from ",x.delta," ");
                              x.delta<-x.delta/2;#cat("to ",x.delta,"\n");
                                 } ## gap too big so decrease x.delta by factor of 5              

              } 
              ## SEVENTH try for x.delta
              if( !(big.ok && small.ok)){
              variate.x<-c(betafixedMode-x.delta,betafixedMode,betafixedMode+x.delta);
              row<-1;
              mymat<-matrix(data=c(variate.x,rep(NA,3)),ncol=2,byrow=TRUE);colnames(mymat)<-c("x","f(x)");## storage for three rows
              for(betafixed in variate.x){## for each of the three initial points mode +/- one other
              gvalue <- .Call("fitabn_marginals",data.df,as.integer(dag.m),as.integer(dim(data.df)[2]),
                            as.integer(var.types),as.integer(max.parents),
                            as.double(mean),as.double(1/sqrt(prec)),as.double(loggam.shape),as.double(1/loggam.inv.scale),
                            as.integer(max.iters),as.double(epsabs),
                            as.integer(verbose),as.integer(error.verbose),
                            grouped.vars,## int.vector of variables which are mixed model nodes
                            group.ids,
                            as.double(epsabs.inner),
                            as.integer(max.iters.inner),
                            as.double(finite.step.size),
                            as.double(hessian.params),
                            ## additional marginal arguments - node, parameter, and relevant modes
                            as.integer(nodeid-1),## C index from 0
                            as.integer(paramid-1),## C index from 0
                            as.double(res.list$modes[[nodeid]]), ## relevant mode estimates
                            as.double(betafixed), ## value to be evaluated e.g. want f(x) given x - betafixed is x
                            as.double(res.list[[names(res.list$modes)[nodeid]]]),## mlik for the node
                            as.integer(max.iters.hessian),as.double(max.hessian.error),
                            as.double(factor.brent),
                            as.integer(maxiters.hessian.brent),
                            as.double(num.intervals.brent),
                  PACKAGE="abn" ## uncomment to load as package not shlib
                  );
                  mymat[row,]<-c(betafixed,gvalue);row<-row+1; 
                  }
                  gvalue<-mymat[,2];## f(x) values
                  big.ok<-small.ok<-FALSE;
                  #cat("diffs:", abs(gvalue[1]-gvalue[2])," ",abs(gvalue[3]-gvalue[2])," 0.02*gvalue=",0.02*gvalue[2]," 0.1*gvalue=",0.1*gvalue[2],"\n");
                  if(   abs(gvalue[1]-gvalue[2])> 0.02*gvalue[2] && abs(gvalue[3]-gvalue[2])> 0.02*gvalue[2] ){big.ok<-TRUE;}
                  if(   abs(gvalue[1]-gvalue[2])< 0.1*gvalue[2] && abs(gvalue[3]-gvalue[2])< 0.1*gvalue[2] ){small.ok<-TRUE;}## ok choices of variate.x
                  if(!big.ok){
                              #cat("increasing 7th time x.delta by *2 from ",x.delta," ");
                              x.delta<-x.delta*2;#cat("to ",x.delta,"\n");
                              }## gap is not big enough so increase x.delta by factor of 5
                   if(!small.ok){
                              #cat("decreasing 7th time x.delta by *2 from ",x.delta," ");
                              x.delta<-x.delta/2;#cat("to ",x.delta,"\n");
                                 } ## gap too big so decrease x.delta by factor of 5              

              }
              ## EIGHTH try for x.delta
              if( !(big.ok && small.ok)){
              variate.x<-c(betafixedMode-x.delta,betafixedMode,betafixedMode+x.delta);
              row<-1;
              mymat<-matrix(data=c(variate.x,rep(NA,3)),ncol=2,byrow=TRUE);colnames(mymat)<-c("x","f(x)");## storage for three rows
              for(betafixed in variate.x){## for each of the three initial points mode +/- one other
              gvalue <- .Call("fitabn_marginals",data.df,as.integer(dag.m),as.integer(dim(data.df)[2]),
                            as.integer(var.types),as.integer(max.parents),
                            as.double(mean),as.double(1/sqrt(prec)),as.double(loggam.shape),as.double(1/loggam.inv.scale),
                            as.integer(max.iters),as.double(epsabs),
                            as.integer(verbose),as.integer(error.verbose),
                            grouped.vars,## int.vector of variables which are mixed model nodes
                            group.ids,
                            as.double(epsabs.inner),
                            as.integer(max.iters.inner),
                            as.double(finite.step.size),
                            as.double(hessian.params),
                            ## additional marginal arguments - node, parameter, and relevant modes
                            as.integer(nodeid-1),## C index from 0
                            as.integer(paramid-1),## C index from 0
                            as.double(res.list$modes[[nodeid]]), ## relevant mode estimates
                            as.double(betafixed), ## value to be evaluated e.g. want f(x) given x - betafixed is x
                            as.double(res.list[[names(res.list$modes)[nodeid]]]),## mlik for the node
                            as.integer(max.iters.hessian),as.double(max.hessian.error),
                            as.double(factor.brent),
                            as.integer(maxiters.hessian.brent),
                            as.double(num.intervals.brent),
                  PACKAGE="abn" ## uncomment to load as package not shlib
                  );
                  mymat[row,]<-c(betafixed,gvalue);row<-row+1; 
                  }
                  gvalue<-mymat[,2];## f(x) values
                  big.ok<-small.ok<-FALSE;
                  #cat("diffs:", abs(gvalue[1]-gvalue[2])," ",abs(gvalue[3]-gvalue[2])," 0.02*gvalue=",0.02*gvalue[2]," 0.1*gvalue=",0.1*gvalue[2],"\n");
                  if(   abs(gvalue[1]-gvalue[2])> 0.02*gvalue[2] && abs(gvalue[3]-gvalue[2])> 0.02*gvalue[2] ){big.ok<-TRUE;}
                  if(   abs(gvalue[1]-gvalue[2])< 0.1*gvalue[2] && abs(gvalue[3]-gvalue[2])< 0.1*gvalue[2] ){small.ok<-TRUE;}## ok choices of variate.x
                  if(!big.ok){
                               #cat("increasing 8th time x.delta by *2 from ",x.delta," ");
                               x.delta<-x.delta*2;#cat("to ",x.delta,"\n");
                              }## gap is not big enough so increase x.delta by factor of 5
                   if(!small.ok){
                               #cat("decreasing 8th time x.delta by *2 from ",x.delta," ");
                               x.delta<-x.delta/2;#cat("to ",x.delta,"\n");
                                } ## gap too big so decrease x.delta by factor of 5              

              }
              ## NINTH try for x.delta
              if( !(big.ok && small.ok)){
              variate.x<-c(betafixedMode-x.delta,betafixedMode,betafixedMode+x.delta);
              row<-1;
              mymat<-matrix(data=c(variate.x,rep(NA,3)),ncol=2,byrow=TRUE);colnames(mymat)<-c("x","f(x)");## storage for three rows
              for(betafixed in variate.x){## for each of the three initial points mode +/- one other
              gvalue <- .Call("fitabn_marginals",data.df,as.integer(dag.m),as.integer(dim(data.df)[2]),
                            as.integer(var.types),as.integer(max.parents),
                            as.double(mean),as.double(1/sqrt(prec)),as.double(loggam.shape),as.double(1/loggam.inv.scale),
                            as.integer(max.iters),as.double(epsabs),
                            as.integer(verbose),as.integer(error.verbose),
                            grouped.vars,## int.vector of variables which are mixed model nodes
                            group.ids,
                            as.double(epsabs.inner),
                            as.integer(max.iters.inner),
                            as.double(finite.step.size),
                            as.double(hessian.params),
                            ## additional marginal arguments - node, parameter, and relevant modes
                            as.integer(nodeid-1),## C index from 0
                            as.integer(paramid-1),## C index from 0
                            as.double(res.list$modes[[nodeid]]), ## relevant mode estimates
                            as.double(betafixed), ## value to be evaluated e.g. want f(x) given x - betafixed is x
                            as.double(res.list[[names(res.list$modes)[nodeid]]]),## mlik for the node
                            as.integer(max.iters.hessian),as.double(max.hessian.error),
                            as.double(factor.brent),
                            as.integer(maxiters.hessian.brent),
                            as.double(num.intervals.brent),
                  PACKAGE="abn" ## uncomment to load as package not shlib
                  );
                  mymat[row,]<-c(betafixed,gvalue);row<-row+1; 
                  }
                  gvalue<-mymat[,2];## f(x) values
                  big.ok<-small.ok<-FALSE;
                  #cat("diffs:", abs(gvalue[1]-gvalue[2])," ",abs(gvalue[3]-gvalue[2])," 0.02*gvalue=",0.02*gvalue[2]," 0.1*gvalue=",0.1*gvalue[2],"\n");
                  if(   abs(gvalue[1]-gvalue[2])> 0.02*gvalue[2] && abs(gvalue[3]-gvalue[2])> 0.02*gvalue[2] ){big.ok<-TRUE;}
                  if(   abs(gvalue[1]-gvalue[2])< 0.1*gvalue[2] && abs(gvalue[3]-gvalue[2])< 0.1*gvalue[2] ){small.ok<-TRUE;}## ok choices of variate.x
                  if(!big.ok){
                               #cat("increasing 9th time x.delta by *2 from ",x.delta," ");
                               x.delta<-x.delta*2;#cat("to ",x.delta,"\n");
                              }## gap is not big enough so increase x.delta by factor of 5
                   if(!small.ok){
                               #cat("decreasing 9th time x.delta by *2 from ",x.delta," ");
                               x.delta<-x.delta/2;#cat("to ",x.delta,"\n");
                                } ## gap too big so decrease x.delta by factor of 5              

              }
 
              ## TENTH try for x.delta
              if( !(big.ok && small.ok)){
              variate.x<-c(betafixedMode-x.delta,betafixedMode,betafixedMode+x.delta);
              row<-1;
              mymat<-matrix(data=c(variate.x,rep(NA,3)),ncol=2,byrow=TRUE);colnames(mymat)<-c("x","f(x)");## storage for three rows
              for(betafixed in variate.x){## for each of the three initial points mode +/- one other
              gvalue <- .Call("fitabn_marginals",data.df,as.integer(dag.m),as.integer(dim(data.df)[2]),
                            as.integer(var.types),as.integer(max.parents),
                            as.double(mean),as.double(1/sqrt(prec)),as.double(loggam.shape),as.double(1/loggam.inv.scale),
                            as.integer(max.iters),as.double(epsabs),
                            as.integer(verbose),as.integer(error.verbose),
                            grouped.vars,## int.vector of variables which are mixed model nodes
                            group.ids,
                            as.double(epsabs.inner),
                            as.integer(max.iters.inner),
                            as.double(finite.step.size),
                            as.double(hessian.params),
                            ## additional marginal arguments - node, parameter, and relevant modes
                            as.integer(nodeid-1),## C index from 0
                            as.integer(paramid-1),## C index from 0
                            as.double(res.list$modes[[nodeid]]), ## relevant mode estimates
                            as.double(betafixed), ## value to be evaluated e.g. want f(x) given x - betafixed is x
                            as.double(res.list[[names(res.list$modes)[nodeid]]]),## mlik for the node
                            as.integer(max.iters.hessian),as.double(max.hessian.error),
                            as.double(factor.brent),
                            as.integer(maxiters.hessian.brent),
                            as.double(num.intervals.brent),
                  PACKAGE="abn" ## uncomment to load as package not shlib
                  );
                  mymat[row,]<-c(betafixed,gvalue);row<-row+1; 
                  }
                  gvalue<-mymat[,2];## f(x) values
                  big.ok<-small.ok<-FALSE;
                  #cat("diffs:", abs(gvalue[1]-gvalue[2])," ",abs(gvalue[3]-gvalue[2])," 0.02*gvalue=",0.02*gvalue[2]," 0.1*gvalue=",0.1*gvalue[2],"\n");
                  if(   abs(gvalue[1]-gvalue[2])> 0.02*gvalue[2] && abs(gvalue[3]-gvalue[2])> 0.02*gvalue[2] ){big.ok<-TRUE;}
                  if(   abs(gvalue[1]-gvalue[2])< 0.1*gvalue[2] && abs(gvalue[3]-gvalue[2])< 0.1*gvalue[2] ){small.ok<-TRUE;}## ok choices of variate.x
                  if(!big.ok){
                              #cat("increasing 10th time x.delta by *2 from ",x.delta," ");
                              x.delta<-x.delta*2;#cat("to ",x.delta,"\n");
                              }## gap is not big enough so increase x.delta by factor of 5
                   if(!small.ok){
                              #cat("decreasing 10th time x.delta by *2 from ",x.delta," ");
                              x.delta<-x.delta/2;#cat("to ",x.delta,"\n");
                                 } ## gap too big so decrease x.delta by factor of 5              

              }
             
              ## BACKUP Use last guess if still not very good
              if( !(big.ok && small.ok)){## initial guesses not good so rebuild variate.x with last guess
              if (verbose) cat("Failed to find a good x-increment so going with last estimate which may not be very good...\n May be better to provide the x-grid manually\n");
              variate.x<-c(betafixedMode-x.delta,betafixedMode,betafixedMode+x.delta);
              row<-1;
              mymat<-matrix(data=c(variate.x,rep(NA,3)),ncol=2,byrow=TRUE);colnames(mymat)<-c("x","f(x)");## storage for three rows
              for(betafixed in variate.x){## for each of the three initial points mode +/- one other
              gvalue <- .Call("fitabn_marginals",data.df,as.integer(dag.m),as.integer(dim(data.df)[2]),
                            as.integer(var.types),as.integer(max.parents),
                            as.double(mean),as.double(1/sqrt(prec)),as.double(loggam.shape),as.double(1/loggam.inv.scale),
                            as.integer(max.iters),as.double(epsabs),
                            as.integer(verbose),as.integer(error.verbose),
                            grouped.vars,## int.vector of variables which are mixed model nodes
                            group.ids,
                            as.double(epsabs.inner),
                            as.integer(max.iters.inner),
                            as.double(finite.step.size),
                            as.double(hessian.params),
                            ## additional marginal arguments - node, parameter, and relevant modes
                            as.integer(nodeid-1),## C index from 0
                            as.integer(paramid-1),## C index from 0
                            as.double(res.list$modes[[nodeid]]), ## relevant mode estimates
                            as.double(betafixed), ## value to be evaluated e.g. want f(x) given x - betafixed is x
                            as.double(res.list[[names(res.list$modes)[nodeid]]]),## mlik for the node
                            as.integer(max.iters.hessian),as.double(max.hessian.error),
                            as.double(factor.brent),
                            as.integer(maxiters.hessian.brent),
                            as.double(num.intervals.brent),
              PACKAGE="abn" ## uncomment to load as package not shlib
              );
              mymat[row,]<-c(betafixed,gvalue);row<-row+1; 
              }}
        ###################################################################################################################
        ### End of successive guesses at initial min. step size
        ## we now have the mode x and its height f(x) and one point either side. 
        ###################################################################################################################
        ## STEP 2. traverse left until we find f(x) less than pdf.min. To do this we use a spline to try and estimate the 
        ##         next x value
        gvalue<-.Machine$double.xmax
         iter<-1;## to avoid taking forever...
         while(gvalue>min.pdf && iter <iter.max){

           gvalue.max<-max(mymat[,2]);## max of f(x)
           x.delta.orig<-x.delta;## backup of stepsize
           max.fact.delta=0.2;## if new x point is more than a factor max.fact.delta (e.g. 0.2) from last evaluated point then stop here
           #cat("evaluating node ",nodeid,": parameter:",paramid," at betafixed=",betafixed," with gvalue=",gvalue,"\n",sep="");  
           ## find the next x value left which differs from the max. gvalue by at least a factor of g.factor, searching in step sizes of
           ##  x.delta subject to the constraint that if we move more than max.fact.delta*last.x then we evaluate here. Avoids big steps. 
           ## this function using spline()
           betafixed<-find.next.left.x(mymat,gvalue.max,g.factor=0.1,x.delta=x.delta,max.fact.delta=max.fact.delta);

           if(have.precision && betafixed<0){## precision terms must not cross zero so try again with smaller x.factor
                                           x.delta<-x.delta/10; if (verbose) cat("Have precision term but jumped to left of zero so skipping...\n");
                                            next; }  
           ## evaluate at next chosen x
           gvalue <- .Call("fitabn_marginals",data.df,as.integer(dag.m),as.integer(dim(data.df)[2]),
                            as.integer(var.types),as.integer(max.parents),
                            as.double(mean),as.double(1/sqrt(prec)),as.double(loggam.shape),as.double(1/loggam.inv.scale),
                            as.integer(max.iters),as.double(epsabs),
                            as.integer(verbose),as.integer(error.verbose),
                            grouped.vars,## int.vector of variables which are mixed model nodes
                            group.ids,
                            as.double(epsabs.inner),
                            as.integer(max.iters.inner),
                            as.double(finite.step.size),
                            as.double(hessian.params),
                            ## additional marginal arguments - node, parameter, and relevant modes
                            as.integer(nodeid-1),## C index from 0
                            as.integer(paramid-1),## C index from 0
                            as.double(res.list$modes[[nodeid]]), ## relevant mode estimates
                            as.double(betafixed), ## value to be evaluated e.g. want f(x) given x - betafixed is x
                            as.double(res.list[[names(res.list$modes)[nodeid]]]),## mlik for the node
                            as.integer(max.iters.hessian),as.double(max.hessian.error),
                            as.double(factor.brent),
                            as.integer(maxiters.hessian.brent),
                            as.double(num.intervals.brent),
              PACKAGE="abn" ## uncomment to load as package not shlib
              );
            if (verbose) cat("left step:",iter," x=",betafixed,#" with min.step.size=",x.delta,
                                                 " f(x)=",gvalue,"\n",sep="");
             mymat<-rbind(mymat,c(betafixed,gvalue));
            iter<-iter+1;
            } ## end of while
            if(iter==1000){ if (verbose) cat("reached max.marg.iter.limit so stopped evaluation and distribution may not be complete!\n");}              
            ## restore x.factor in case changed during precision evaluation.
            x.delta<-x.delta.orig;
        #####################################################################################################################
        ## STEP 3. as above but traverse right 
        ## 
        iter<-1;## to avoid taking forever...
        gvalue<-.Machine$double.xmax
        if (verbose) cat("gvalue=",gvalue," min.pdf=",min.pdf," iter.max=",iter.max,"\n");
         while(gvalue>min.pdf && iter <iter.max){

           gvalue.max<-max(mymat[,2]);## max of f(x)
           x.delta.orig<-x.delta;## backup of stepsize
           max.fact.delta=0.2;## if new x point is more than a factor max.fact.delta (e.g. 0.2) from last evaluated point then stop here
           #cat("evaluating node ",nodeid,": parameter:",paramid," at betafixed=",betafixed," with gvalue=",gvalue,"\n",sep="");  
           ## find the next x value left which differs from the max. gvalue by at least a factor of g.factor, searching in step sizes of
           ##  x.delta subject to the constraint that if we move more than max.fact.delta*last.x then we evaluate here. Avoids big steps. 
           ## this function using spline()
           betafixed<-find.next.right.x(mymat,gvalue.max,g.factor=0.1,x.delta=x.delta,max.fact.delta=max.fact.delta);
 
           ## evaluate at next chosen x
           gvalue <- .Call("fitabn_marginals",data.df,as.integer(dag.m),as.integer(dim(data.df)[2]),
                            as.integer(var.types),as.integer(max.parents),
                            as.double(mean),as.double(1/sqrt(prec)),as.double(loggam.shape),as.double(1/loggam.inv.scale),
                            as.integer(max.iters),as.double(epsabs),
                            as.integer(verbose),as.integer(error.verbose),
                            grouped.vars,## int.vector of variables which are mixed model nodes
                            group.ids,
                            as.double(epsabs.inner),
                            as.integer(max.iters.inner),
                            as.double(finite.step.size),
                            as.double(hessian.params),
                            ## additional marginal arguments - node, parameter, and relevant modes
                            as.integer(nodeid-1),## C index from 0
                            as.integer(paramid-1),## C index from 0
                            as.double(res.list$modes[[nodeid]]), ## relevant mode estimates
                            as.double(betafixed), ## value to be evaluated e.g. want f(x) given x - betafixed is x
                            as.double(res.list[[names(res.list$modes)[nodeid]]]),## mlik for the node
                            as.integer(max.iters.hessian),as.double(max.hessian.error),
                            as.double(factor.brent),
                            as.integer(maxiters.hessian.brent),
                            as.double(num.intervals.brent),
              PACKAGE="abn" ## uncomment to load as package not shlib
              );
            if (verbose) cat("right step:",iter," x=",betafixed,#" ,min.step.size=",x.delta,
                                                   " f(x)=",gvalue,"\n",sep="");
             mymat<-rbind(mymat,c(betafixed,gvalue));
            iter<-iter+1;
            } ## end of while
            if(iter==1000){if (verbose) cat("reached max.marg.iter.limit so stopped evaluation and distribution may not be complete!\n");}              
 
            ## final step before saving data into list is to make sure it is ordered according to increasing x
            mymat<-mymat[order(mymat[,1]),];
            node.marginals[[curnom]]<-mymat;

             } ## parameter

      
       
       marginals[[index]]<-node.marginals;index<-index+1; }  ## check whether not already computed in INLA
     
     } ## node
    
     } ## end of if single node or all

     #cat("f(",betafixed,")=",marg.res,"\n");   
if( is.null(marginal.node)) {## not in single node case   
    names(marginals)<-colnames(dag.m)[which(INLA.marginals==FALSE)]
    res.list[["marginals"]]<-marginals;
} else {## single parameter manual case
        res.list[["marginals"]]<-marginals;
}

return(res.list);


}

#############################################################################
## attempt to find the next x evaluation point using spline extrapolation
## traversing left from mode
#############################################################################
 find.next.left.x<-function(mat.xy,g.max,g.factor,x.delta,max.fact.delta){
 
   xleft<-min(mat.xy[,"x"]);## current left hand endpoint  
   orig.xleft<-xleft;
   do<-TRUE;       
   while(do){
     ## keep decreasing xleft until f(x) reduces to next point of evaluation
    x.out=xleft-abs(x.delta);## size of vertical grid to the left is abs(x.delta)
    if(abs(x.out-orig.xleft)>max.fact.delta*x.delta){do<-FALSE;} ## if current value of x is a distance further away than original value then stop here
    next.try<-spline(mat.xy[,"x"],mat.xy[,"f(x)"],xout=x.out);## get value of f(x) at x.factor of current xleft 
    if(abs(next.try$y-g.max)>g.factor*g.max){## if current value of x gives a value for g=f(x) which is at least a min vertical distance from current then finish
    return(next.try$x);}

    ## get to here need to keep moving left at the estimated value of g is not sufficiently different from the current one
    xleft<-next.try$x;## try a lower xleft 
            } ## end of while
 
   return(xleft);
  }
#############################################################################
## attempt to find the next x evaluation point using spline extrapolation
## traversing right from mode
#############################################################################
 find.next.right.x<-function(mat.xy,g.max,g.factor,x.delta,max.fact.delta){
 
   xright<-max(mat.xy[,"x"]);## current left hand endpoint  
   orig.xright<-xright;
   do<-TRUE;       
   while(do){
     ## keep increasing xright until f(x) reduces to next point of evaluation
    x.out=xright+abs(x.delta);## size of vertical grid to the left is abs(x.delta)
    if(abs(x.out-orig.xright)>max.fact.delta*x.delta){do<-FALSE;} ## if current value of x is a distance further away than original value then stop here
    next.try<-spline(mat.xy[,"x"],mat.xy[,"f(x)"],xout=x.out);## get value of f(x) at x.out
    if(abs(next.try$y-g.max)>g.factor*g.max){## if current value of x gives a value for g=f(x) which is at least a min vertical distance from current then finish
    return(next.try$x);}

    ## get to here need to keep moving right at the estimated value of g is not sufficiently different from the current one
    xright<-next.try$x;## try a higher xright
            } ## end of while
 
   return(xright);
  }
