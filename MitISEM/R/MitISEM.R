# Function to approximate a 'kernel' using an adaptive mixture of 
# multivariate t densities by the IS weighted EM algorithm
# 
# inputs:
#   KERNEL   : [function] which computes the kernel. Must be vectorized
#   mu0      : [kx1 vector] initial value for the location of the first component
#   Sigma0   : [kxk matrix] scaling matrix of the first component (default: NULL, estimated by MitISEM)
#   df0      : [double>0] degrees of freedom of the first component (default: 1)
#   mit0     : [list] with initial mixture component, default: NULL, optimized using mu0,Sigma0,df0
#   control  : [list] control parameters with the following components:
#     N           : [integer>100] number of draws used in the simulation (default: 1e5)
#     robust.N    : [logical] if 'TRUE' get 'N' draws with finite KERNEL values (default: TRUE)
#     Hmax        : [integer>0] maximum number of components (default: 10)
#     StopMethod  : [string], 'CV' or 'AR' defines type of stopping criteria for MitISEM approximation, see 'details' 
#     CVtol       : [double in (0,1)] convergence criteria for CV (default: 0.1, i.e. 10%)
#                    used only if 'StopMethod=='CV''
#     ARtol      : [double in (0,1)] convergence criteria for Mixture probabilities (default: 0.1, i.e. 10%)
#                    used only if 'StopMethod=='AR''
#     trace       : [logical] output printed during the fitting (default: FALSE)
#     trace.init  : [logical] output printed of the optimizer (default: 0, i.e. no output)
#     maxit.init  : [double] maximum number of iterations in the optimization (default: 1e4)
#     reltol.init : [double] relative tolerance in the optimization (default: 1e-8)
#     maxit.EM    : [integer>0] max. number of iterations for the EM algorithm, default: 1000
#     tol.EM      : [double>0] tolerance for EM steps' convergence, default: 0.001
#     trace.EM    : [logical] tracing information on the fitting procedure of IS_EM
#                   (Default: FALSE; i.e., no tracing information)
#     optim.df    : [logical] default: TRUE, df are optimized, df fixed at initial level otherwise
#                   (following inputs active only if df are optimized: optim.df=TRUE)
#     inter.df    : [2-vector] range of search values for df optimization. Default: c(0.01,30) (min and max degrees of freedom)
#     tol.df      : [double>0] tolerance for degree of freedom optimization step (default 0.0001220703)
#     maxit.df    : [integer>0] max. # of iterations for df optimization (default 1000)
#     trace.df    : [logical] tracing information on the fitting procedure of degrees of freedom. 
#                   (Default: FALSE, i.e., no tracing information.
#     tol.pr      : [double in [0,1)] minimum probability required to keep mix. components
#                     (default: 0, component probability does not alter mixture)
#     ISpc        : [double in (0,1)] fraction of draws to construct new component (default: 0.1)
#     Pnc         : [double in (0,1)] initial probability of the new component (default: 0.1)
#   ...         : additional parameters used by 'KERNEL'
# outputs:
#   mit        : [list] with the following components:
#      p       : [Hx1 vector] of probabilities
#      mu      : [Hxk matrix] of location vectors (in rows)
#      Sigma   : [Hxk^2 matrix] of scale matrices (in rows)
#      df      : [Hx1 vector] of degrees of freedom
#   CV         : [Hx1 vector] of coefficient of variation
#   time       : [double] processed time
#   summary    : [data.frame] information on construction of components, time and CV
#
# author : Nalan Basturk
# date   : 20120912

MitISEM <- function(KERNEL,mu0,Sigma0=NULL,df0=1,mit0=NULL,control=list(),...){
   # Check inputs         
   if(missing(KERNEL)) 
      stop ("'KERNEL' is missing in 'MitISEM'")
   KERNEL <- match.fun(KERNEL)
   if(!any(names(formals(KERNEL))=="log"))
      stop ("'KERNEL' MUST have the logical argument 'log' returning log(KERNEL) if 'log=TRUE'")  
   if(missing(mu0) & is.null(mit0))
      stop ("'mu0' or 'mit0' has to be supplied in 'MitISEM'")
   if(is.null(mit0)){
      if(missing(mu0) | !is.vector(mu0))
        stop("if 'mit0' is not supplied, 'mu0' must be a vector of size 'k'")
      if(!is.null(Sigma0) & !is.matrix(Sigma0))
        stop("'Sigma0' must be a 'kxk' matrix")
      if(!is.null(Sigma0) & any(dim(Sigma0)!=length(mu0)))
        stop("Dimensions of 'Sigma0' and 'mu0' do not match")
   }
   if(df0 <= 0){
      df0=1; warning("'df0' must be above 0, default value '1' is used instead")
   }
   if(!is.null(mit0)){
      if(!isMit(mit0))
       stop ("initial mixture 'mit0' not correctly specified in MitISEM")
      if(any(mit0$df<=0))
       stop ("initial mixture 'mit0' must have positive degrees of freedom")
   }
   
   # Set/check control parameters
   con <- list(N=1e4,robust.N=TRUE,Hmax=10,StopMethod='CV', CVtol=0.1, ARtol = 0.1,
            trace=FALSE, trace.init=FALSE, maxit.init=1e4, reltol.init=.Machine$double.eps^0.25,
            maxit.EM = 1e3, tol.EM = 0.001, trace.EM=FALSE,
            optim.df=TRUE, inter.df=c(0.01,30), tol.pr=0, tol.df=.Machine$double.eps^0.25,
            maxit.df=1e3, trace.df=TRUE,ISpc = 0.1, Pnc=0.1)
   con[names(control)] <- control
   control.ISEM <- list(maxit.EM = con$maxit.EM, tol.EM = con$tol.EM,trace.EM=con$trace.EM,
                    optim.df=con$optim.df,inter.df=con$inter.df,
					tol.pr = con$tol.pr,
                    tol.df=con$tol.df,maxit.df=con$maxit.df,trace.df=FALSE)
   control.init <- list(trace=con$trace.init,maxit=con$maxit.init,reltol=con$reltol.init)
   if(con$N<100)
      stop ("'N', number of parameter draws is too small.")
   if(con$Hmax<1)
      stop ("'Hmax' must be at least '1'")
   if(!is.logical(con$optim.df))
      stop ("'optim.df' has to be logical")
   if(con$CVtol<=0 | con$CVtol>=1)
      stop ("'CVtol' must lie in (0,1)")
   if(con$ARtol<=0 | con$ARtol>=1)
      stop ("'CVtol' must lie in (0,1)")     
   if(con$tol.EM <= 0)
      stop("tol.EM must be double > 0");
   if(!is.logical(con$trace.EM)) 
      stop("'trace.EM' must be logical");
   if(!is.logical(con$optim.df)) 
      stop("'optim.df' must be logical");
   if(!all(range(con$inter.df)==con$inter.df) | length(con$inter.df)!=2 | any(con$inter.df <= 0))
      stop("'iter.df' must be vector of length 2, with positive, increasing values");
   if(con$tol.pr < 0 | con$tol.pr >= 1)
      stop("tol.pr must be double in [0,1)");
   if(con$tol.df <= 0)
      stop("tol.df must be double > 0");
   if(!is.logical(con$trace.df)) 
      stop("'trace.df' must be logical");
   if(con$ISpc<=0 |con$ISpc>=1)
      stop ("'ISpc' must be a double in (0,1)")
   if(con$Pnc<=0 |con$Pnc>=1)
      stop ("'Pnc' must be a double in (0,1)")
   if(!is.logical(con$trace))
      stop ("'trace' must be logical")
   if(!is.logical(con$robust.N))
      stop ("'robust.N=TRUE' must be logical")     
     
   N      <- con$N
   Hmax   <- con$Hmax
   CVtol  <- con$CVtol
   ARtol  <- con$ARtol
   m_stop <- con$StopMethod
   trace  <- con$trace
   
   # STEP 1: get/define initial mit density shape, scale and degree of freedom
   if(is.null(mit0)){
     mit.init <- list(p=NULL, mu=NULL, Sigma=NULL, df=df0)
	 control.init$hessian=is.null(Sigma0)
     if (is.null(Sigma0)) { # optimize initial scale(s)
        initopt <- fn.initoptim(KERNEL=KERNEL, mu0, control=control.init,...)
     }else {  # user defined sscale
        initopt <- list(mu=mu0, Sigma=Sigma0, method="USER", time=0)
     }
     mit.init$p     <- 1  
     mit.init$mu    <- matrix(initopt$mu, nrow=1)
     mit.init$Sigma <- matrix(initopt$Sigma, nrow = 1)
   }else{
     mit.init = mit0 # user defined density
     initopt <- list(method="USER", time=0)
   }
   # get draws and IS weights from initial mixture
   tmp      <- fn.rmvgt_robust(robustify=con$robust.N,N,mit=mit.init,KERNEL,log=TRUE,...)
   theta    <- tmp$theta 
   lnk      <- tmp$lnk
   lnd      <- dmvgt(theta, mit.init, log=TRUE)
   w        <- fn.ISwgts(lnk, lnd)
   
   # update scale and location using IS-EM, fixing 'df' 
   cont.init <- control.ISEM
   cont.init$optim.df = FALSE
   mit.init <- fn.optimt(theta,mit.init,w,control=cont.init)$mit
   tmp      <- fn.rmvgt_robust(robustify=con$robust.N,N,mit=mit.init,KERNEL,log=TRUE,...)
   theta    <- tmp$theta 
   lnk      <- tmp$lnk
   lnd      <- dmvgt(theta, mit.init, log=TRUE)
   w        <- fn.ISwgts(lnk, lnd)
   H        <- length(mit.init$p) # number of components (1 if mit0=NULL)
   
   # stopping criteria
   stopcr   <- fn.stopMit(method=m_stop,w,CV_last=100,CVtol,AR_last=100,ARtol)
   CV       <- stopcr[1]
   AR       <- stopcr[2]
  
   # summary
   nsummary <- c("H","METHOD","TIME","CV") 
   summary  <- data.frame(H, initopt$method, initopt$time,  CV)
   names(summary) = nsummary
   if (trace)  
     print(summary,row.names=FALSE);
   # STEP 2: optimize mixture using IS weighted EM and get draws from the new mit
   # optimize mode scale and df (df can be fixed by user)
   ptm        <- proc.time()[3]
   optimt     <- fn.optimt(theta,mit.init,w,control=control.ISEM)
   mit.new    <- optimt$mit
   # get draws and log kernel evaluation 
   tmp        <- fn.rmvgt_robust(robustify=con$robust.N,N,mit=mit.new,KERNEL,log=TRUE,...)
   theta.new  <- tmp$theta 
   lnk        <- tmp$lnk
   lnd        <- dmvgt(theta.new, mit.new, log=TRUE)
   w          <- fn.ISwgts(lnk, lnd)
   stopcr     <- fn.stopMit(method=m_stop,w,CV_last=100,
                          CVtol,AR_last=100,ARtol)
   CV         <- c(CV,stopcr[1]) # initial coefficient of variation
   AR         <- c(AR,stopcr[2]) # initial acceptance rate
   summary.n  <- data.frame(H, "IS-EM", as.numeric(proc.time()[3]-ptm), stopcr[1])
   names(summary.n)=nsummary
   summary = rbind(summary,summary.n)
   if (trace)  print(summary,row.names=FALSE);
   
   # STEP 3: add more mixture components in 'mit' until convergence
   hstop <- FALSE;
   while (H<Hmax & hstop==FALSE){
      ptm <- proc.time()[3]
      H    = H + 1
      # select largest IS weights and corresponding draws
      ind.w    <- fn.select(w,con$ISpc)
      theta.nc <- theta.new[ind.w,]
      w.nc     <- w[ind.w]

      # compute new component's mode and scale from IS weights
      mit.nc       <- list(p=con$Pnc,df=1)
      if(!control.ISEM$optim.df)
        mit.nc$df = df0
      tmp          <- fn.optIS(theta.nc,w.nc)
      mit.nc$mu    = tmp$mu
      mit.nc$Sigma = tmp$Sigma

      # combine old 'mit' and new t component
      mit.new   <- fn.update.mit(mit.new,mit.nc)
      # get draws and log kernel evaluation from new 'mit'
      tmp       <- fn.rmvgt_robust(robustify=con$robust.N,N,mit=mit.new,KERNEL,log=TRUE,...)
      theta.new <- tmp$theta 
      lnk       <- tmp$lnk
      lnd       <- dmvgt(theta.new, mit.new, log=TRUE)
      w         <- fn.ISwgts(lnk, lnd)
      
      # update mode/scale/df of all mixture components
      mit.new <- fn.optimt(theta.new,mit.new,w,control.ISEM)$mit
      H        = length(mit.new$p)

      # draw new parameters from 'mit', evaluate new IS weights and convergence
      tmp       <- fn.rmvgt_robust(robustify=con$robust.N,N,mit=mit.new,KERNEL,log=TRUE,...)
      theta.new <- tmp$theta 
      lnk       <- tmp$lnk   
      lnd       <- dmvgt(theta.new, mit.new, log=TRUE)
      w         <- fn.ISwgts(lnk, lnd)
      stopcr    <- fn.stopMit(method=m_stop,w,CV_last=CV[length(CV)],
                              CVtol,AR_last=AR[length(CV)],ARtol)
      CV        <- c(CV,stopcr[1])
      AR        <- c(AR,stopcr[2])
      
      if(H > 1)
        hstop  <- as.logical(stopcr[3]) 
       
      summary.n <- data.frame(H, "IS-EM", as.numeric(proc.time()[3]-ptm),  CV[length(CV)])
      names(summary.n)=nsummary
      summary <- rbind(summary, summary.n)
      row.names(summary) <- NULL
      if (trace)  
	    print(summary.n,row.names=FALSE);
    }
    return(list(mit=mit.new,CV=CV,time=sum(summary$TIME),summary=summary))
}
# 'fn.optIS' updates mode and scale of the multivariate Student-t matrix from IS weights
# inputs:
#    theta: (T* x k) matrix of draws with highest IS weights
#    w    : vector of size T* of IS weights corresponding to 'theta'
# outputs: list containing
#    mu   : optimized mode of the student t component
#    Sigma: optimized scale of the student t component
fn.optIS <- function(theta,w){
    theta = as.matrix(theta)
    w     = exp(log(w) - log(sum(w))) # normalized weights    
    mu    <- apply(theta,2,function(x)(sum(x*w)))
    er    <- t(apply(theta,1,function(x)(x-mu))) 
    wer   <- w * er
    Sigma <- t(er) %*% wer    
    mu    <- matrix(mu,nrow=1)     #return mode and scale in MitISEM format
    Sigma <- matrix(Sigma,nrow=1)    
    return(list(mu=mu,Sigma=Sigma)) 
}
# 'fn.update.mit' COMBINES OLD MIXTURE AND THE NEW COMPONENT
# inputs:
#   mit     : last mit density of H student-t components
#   mit.nc  : new mit density of 1 student-t component #   outputs:
# outputs:
#   mit.m   : mixture of t (mit) density combining all components
fn.update.mit <- function(mit,mit.nc){
  mit.m       <- list(p=NULL,mu=NULL,Sigma=NULL,df=NULL)
  Pnc         <- mit.nc$p           #probability of the new component
  mit.m$p     = c((1-Pnc)*mit$p, Pnc)
  mit.m$mu    = rbind(mit$mu,mit.nc$mu)
  mit.m$Sigma = rbind(mit$Sigma,mit.nc$Sigma)
  mit.m$df    = c(mit$df,mit.nc$df)
  (mit.m)   
}
# 'fn.select' gets index of draws with largest IS weights
# inputs:
#   w      : [Nx1] vector of IS weights
#   pc     : [double in (0,1)] percentage for highest weighted IS draws to use
# outputs:
#   ind    : [Nx1] vector of index of highest IS weight draws
fn.select<-function(w,pc){
   N    <- length(w)
   xr   <- round(N*pc)-1
   ind  <- order(w)[(N-xr):N]
   (as.vector(ind))
}

