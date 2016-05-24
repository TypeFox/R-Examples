# This function optimizes the location and the scale of the initial 
# multivariate (k-variate) t distribution in MitISEM
#
# inputs:
#   KERNEL   : [function] to be approzimated
#   mu0      : [vector size k] starting value for location
#   control  : [list] of optimization parameters (see optim)
#   ...      : additional inputs used by 'KERNEL'
# outputs: [list] with the following components:
#   mu       : [vector size k] of location parameters
#   Sigma    : [vector size k^2] of scale matrices (in vector form)
#   method   : [string] method used in the optimization: "BFGS" or "Nelder-Mead"
#   time     : [double] time
#
# note: Scale variable is optimized if not user-defined
# author : Nalan Basturk
# date   : 20120912
        
fn.initoptim <- function(KERNEL, mu0, control, ...)
{
    hessian = control$hessian
    control<-list(trace=control$trace,maxit=control$maxit,
            reltol=control$reltol,fnscale=-1)
    ptm <- proc.time()[3]
    k <- length(mu0)
    Sigma <- NULL
        # try two optimziation methods
        method <- "BFGS"
        opt <- optim(par=mu0, fn=KERNEL, method=method,control=control,
                 hessian=hessian, log=TRUE,...)
        if(opt$convergence > 1)
          method <- "Nelder-Mead"
          opt <- optim(par=mu0, fn=KERNEL, method=method,control=control,
                 hessian=hessian, log=TRUE,...)
        # if optimized, check if scale is pds
        if(hessian){
          Sigma <- -solve(opt$hessian)
          r     <- fn.isPDS(Sigma)
          if(r) 
            stop("problem in initial optimization, try different initial values")
        }
        list(mu=as.vector(opt$par), Sigma=Sigma, method=method,time=(proc.time()[3]-ptm))
}
