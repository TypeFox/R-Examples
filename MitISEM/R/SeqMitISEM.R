# Function to approximate a 'kernel' using adaptive mixture of 
# multivariate t densities optimized by the IS weighted EM algorithm
# at increasing data samles
#
# inputs:	
#   data             : [Txm matrix or T vector] of data. 
#                      All data to pass to KERNEL must be provided in this argument  
#   KERNEL           : [function] which computes the kernel. Must be vectorized.
#                      KERNEL should have argument 'data' (Txm) matrix for T observations and m data series
#                      KERNEL should have argument 'log', logical indicating log density to return
#   Sigma0           : [kxk matrix] scaling matrix of the first component 
#                      (default: NULL, estimated by MitISEM)
#   df0              : [double>0] degrees of freedom of the first component (default: 1)
#   control.MitISEM  : [list] control parameters to pass to 'MitISEM', see 'MitISEM' for details and defaults
#   control.seq      : [list] containing sequential optimization controls:
#      T0            : [integer] number of observations to use in initial MitISEM, default: 0.5
#      tau           : [l vector] observations to include in the sample, 
#                      default: 1, a single observation added to initial sample
#                      T0 + max(tau)) must be less than the data size
#      trace         : [logical] indicator to trace approximations
#   ...              : additional parameters used by 'KERNEL'
# outputs:
#   mit.init         : [list] MitISEM approximation for the initial data sample, see 'isMit'
#   CV.init          : [double] CoV from first approximation
#   CV.SeqMitISEM    : [l vector] CoV from each sequential approximation
#   mit.pred         : [list length l] MitISEM approximation from each sequential approximation
#                      each component is a Mit density
#   summary          : [matrix] summary of approximations
#
# author : Nalan Basturk
# date   : 20120912

SeqMitISEM <- function(data,KERNEL,mu0,Sigma0=NULL,df0=1,control.MitISEM=list(),control.seq=list(),...){
   if(missing(data))
      stop ("data 'data' is missing in 'SeqMitISEM'")
   if(!is.vector(data) & !is.matrix(data))
      stop ("data 'data' must be vector or matrix in 'SeqMitISEM'")
   if(missing(KERNEL))
      stop ("'KERNEL' is missing in 'SeqMitISEM'")
   KERNEL <- match.fun(KERNEL)
   if(!any(names(formals(KERNEL))=="data"))
     stop ("'KERNEL' must have the argument 'data'")  
   if(!any(names(formals(KERNEL))=="log"))
     stop ("'KERNEL' must have the argument 'log'")  
   if(missing(mu0))
      stop ("'mu0' is missing in 'SeqMitISEM'")
   if(is.null(control.MitISEM$N))
      control.MitISEM$N = 1e4
   if(is.null(control.MitISEM$Hmax))
      control.MitISEM$Hmax = 10
   if(is.null(control.MitISEM$robust.N))
      control.MitISEM$robust.N = TRUE
   if(is.null(control.MitISEM$CVtol))
      control.MitISEM$CVtol = 0.1
   if(is.null(control.MitISEM$ARtol))
      control.MitISEM$ARtol = 0.1     

   N = control.MitISEM$N # number IS draws
   if(is.vector(data))
     data <- as.matrix(data)
   Tend    <- nrow(data)
   con.seq <- list(T0=round(Tend/2),tau=1,tol.seq=0.2,metod=1,trace=TRUE) # control arguments 
   con.seq[names(control.seq)] <- control.seq

   if(max(con.seq$T0+max(con.seq$tau)) > Tend)
      stop("data partitions are not well-defined in SeqMitISEM")
  
   # MitISEM approximation for the kernel evaluated at initial sample
   data_train <- as.matrix(data[1:con.seq$T0,])  
   MitISEM0   <- MitISEM(KERNEL=KERNEL,mu0=mu0,Sigma0=Sigma0,df0=df0,
                      control=control.MitISEM,data=data_train,...)
   H0          <- length(MitISEM0$mit$p)
   CV0         <- MitISEM0$CV[H0+1] # CV in the 2nd step of SeqMitISEM
   MitISEM0    <- list(mit=MitISEM0$mit,H=H0,CV=CV0)

   # Sequential MitISEM for the kernel evaluated at each 'extended sample'
   summary       <- NULL
   mit.pred      <- list()
   CV.SeqMitISEM <- c()
   for (Tstep in con.seq$tau){
     # Select training and estimation sample
     ind      <- Tstep-min(con.seq$tau)+1 # storage index
     Tlast    <- con.seq$T0+Tstep
     data_new <- as.matrix(data[1:Tlast,])
     
     # Apply Sequential MitISEM steps 
     out <- fn_sub.SeqMitISEM(data_new,KERNEL,MitISEM0,
                              control.MitISEM,con.seq$tol.seq,...)
     mit.pred[[ind]]    =  out$mit  
     CV.SeqMitISEM[ind] <- out$CV
     
     # Update initial mit density and COV if step3 is reached
     if(out$update)
       MitISEM0 <- list(mit=out$mit, H=length(out$mit$p),CV=out$CV)      
     
     # summary of results
     summary = rbind(summary,out$summary)
     if(con.seq$trace)
       print(out$summary);
     tmp <- list(mit.init=MitISEM0,CV.SeqMitISEM=CV.SeqMitISEM,
                 mit.pred=mit.pred,summary=summary)
   }
   return(list(mit.init=MitISEM0,CV.init=CV0,CV.SeqMitISEM=CV.SeqMitISEM,
           mit.pred=mit.pred,summary=summary))
}

fn_sub.SeqMitISEM <- function(data_new,KERNEL,MitISEM0,control.MitISEM,tol,...){
  N          <- control.MitISEM$N              # number of IS draws
  Hmax       <- control.MitISEM$Hmax           # max number of mix. components
  nadapt     <- nextend <- nreuse <- 0         # number of adapted/reused/extended mixtures
  mit0       <- MitISEM0$mit                   # initial MitISEM candidate (from previous sample)
  H0         <- MitISEM0$H                     # initial number of mixtures
  CV0        <- MitISEM0$CV                    # coeff. of variation from initial sample
  CV.reuse   <- CV.ad <- CV.ex  <- "-"         # coeff. of variation from adap/reusing/extending
  
  # draws and density from initial MitISEM candodate using new data sample
  tmp        <- fn.rmvgt_robust(robustify=control.MitISEM$robust.N,N,
                                mit=mit0,KERNEL,log=TRUE,data=data_new,...)
  theta.new  <- tmp$theta 
  lnk        <- tmp$lnk
  lnd        <- dmvgt(theta.new, mit0, log=TRUE)
  w          <- fn.ISwgts(lnk, lnd)
 
  stopcr     <- fn.stopMit(method='CV',w,CV_last=100,
                           control.MitISEM$CVtol,AR_last=100,control.MitISEM$ARtol)
  CV.reuse   <- stopcr[1] # coefficient of variation 'reusing' the candidate
  
  # comparison of CV for the training sample and the new sample
  reuse <- (abs((CV.reuse-CV0)/CV0) <= tol)  # logial to 'reuse' previous candidate
  mit.new = mit0
  
  CV.opt <- mit.opt <- NULL;
  if(reuse){
    nreuse = nreuse+1
    CV.opt = CV.reuse
    mit.opt = mit0
  }else{
    # adapted mixture density only with an IS EM step
    control.ad       <- control.MitISEM 
    control.ad$Hmax  = H0 # do not increase the number of components
    MitISEM.ad       <- MitISEM(KERNEL=KERNEL,mit0=mit0,control=control.ad,data=data_new,...)
    H.ad   <- length(MitISEM.ad$mit$p)
    mit.ad <- MitISEM.ad$mit
    CV.ad  <- MitISEM.ad$CV[length(MitISEM.ad$CV)] 
    H0     = H.ad
    # CV comparison using the adapted density
    adapt <- (abs((CV.ad-CV.reuse)/CV.reuse) <= tol) # logial to 'adapt' previous candidate
    if(adapt | H0==Hmax ){
      nadapt=nadapt+1
      CV.opt = CV.ad
      mit.opt = mit.ad
    }else{
      # extending the candidate given the 'adapted' candidate as starting density 
      MitISEM.ex <- MitISEM(KERNEL=KERNEL,mit0=mit.ad,control=control.MitISEM,data=data_new,...)
      mit.ex <- MitISEM.ex$mit
      CV.ex  <- MitISEM.ex$CV[length(CV.ex)]
      nextend = nextend+1
      CV.opt = CV.ex
      mit.opt = mit.ex
    }
  }
  H.opt = length(mit.opt$p)
  
  # summary of the sequential step
  summary = data.frame(H0,H.opt,CV.opt,nreuse,nadapt,nextend)
  nsummary = c("#t in training","#t in optimal","CV","# reuse","# adapt","# extend")
  names(summary) = nsummary
  
  return(list(mit=mit.opt,CV=CV.opt,H=H.opt,update=(!reuse),summary=summary))
}
