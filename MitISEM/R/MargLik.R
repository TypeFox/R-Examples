# Function to evaluate marginal likelihood from Importance Sampling
#
# inputs:
#     N       : [integer > 100] number of draws for IS, default N = 1e4
#     mit     : [list] containing candidate mixture of t density, see 'isMit'
#     KERNEL  : [function] calculating the exact posterior density. KERNEL must have the argument:
#                'log', indicator for returning log-density
#     ...     : extra arguments to pass to 'KERNEL'
# outputs: 
#     ML.mean : [double] estimated marginal likelihood * 10^scale
#                (scaling may be necessary for numerical accuracy, typically with very low ML values)
#     ML.NSE  : [double] Numerical Standard Error from Importance Sampling
#     scale   : [integer > 0] providing scaling of the ML
#
# note: 'scale' is added since ML values can be very small. 
#       'KERNEL' has to be the exact posterior density
#
# author : Nalan Basturk
# date   : 20120912

MargLik <- function(N=1e4,mit,KERNEL,...){
  # check input
  if(missing(mit))
    stop("''mit' must be defined in 'MargLik'")
  if(N <= 100)
    stop("Number of draws 'N' must be above 100 in 'MargLik'")
  if(!isMit(mit))
    stop("''mit' must be lists for mixture of t densities in 'MargLik', see 'isMit'")
  if(missing(KERNEL))
    stop ("'KERNEL' is missing in 'MargLik'")
  KERNEL <- match.fun(KERNEL)
  if(!any(names(formals(KERNEL))=="log"))
    stop ("'KERNEL' must have the argument 'log' in 'MargLik'")  
  
  # draws from candidate density
  theta <- as.matrix(rmvgt(N,mit))
  # log target and candidate densities given draws
  lnk <- apply(theta,1,FUN=KERNEL,log=TRUE,...)
  lnq <- apply(theta,1,FUN=dmvgt,mit=mit,log=TRUE)
    
  # errors if target evaluations are weird
  if(all(lnk==-Inf) | all(lnk==lnk[1]))
    stop("Target evaluations are infinite or constant try more draws 'N' for IS")
    
  # robustify/scale marginal likelihood if weights are too small
  scale = 0; w=0
  while (floor(mean(w))==0){
    lnw = lnk - lnq + scale*log(10) 
    w   = exp(lnw)
    scale = scale + 1
  }
  scale = scale - 1
  
  # ML and NSE from IS weights
  ML.mean = mean(w)
  ML.NSE  = sqrt(var(w)/N)
  return(list(ML.mean=ML.mean,ML.NSE=ML.NSE,scale=scale))
}  
