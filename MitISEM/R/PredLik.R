# Function to calculate predictive likelihoods from Importance Sampling
# using atraining sample and full sample of observations
# and two MitISEM candidate densities for IS
#
# inputs:
#     N       : [integer > 100] number of draws for IS
#     mit.fs  : [list] containing candidate mixture of t density for the 'full sample', see 'isMit'
#     mit.ss  : [list] containing candidate mixture of t density for the 'subsample', see 'isMit'
#     KERNEL  : [function] calculating target density. Should have arguments
#                data (vector or matrix data input)
#     data.fs : [T1xk matrix or vector of size T1] full data
#     data.ss : [T2xk matrix or vector of size T2] full data
#               Note: T2 < T1 should hold
#     ...     : extra arguments to pass to 'KERNEL'
# outputs: 
#     PL      : [double] predictive likelihood x 10^scale
#             (calculated from marginal likelihoods (full sample and subsample))
#                (scaling may be necessary for numerical accuracy, typically with very low ML values)
#     scale   : [integer > 0] providing scaling of the PL
#
# note: 'KERNEL' can be proportional to the posterior density
#       'scale' is added since Predictive Likelihood values can be very small. 
# author : Nalan Basturk
# date   : 20120912

PredLik <- function(N=1e4,mit.fs,mit.ss,KERNEL,data.fs,data.ss,...){
  # check input
  if(N <= 100)
    stop("Number of draws 'N' must be above 100 in 'PredLik'")
  if(missing(mit.fs) | missing(mit.ss))
    stop("''mit.fs' and 'mit.ss' must be defined in 'PredLik'")
  if(missing(data.fs) | missing(data.ss))
    stop("''data.fs' and 'data.ss' must be defined in 'PredLik'")
  if(!isMit(mit.fs) | !isMit(mit.ss))
    stop("''mit.fs' and/or 'mit.ss' must be lists for mixture of t densities in 'PredLik', see 'isMit'")
  if(missing(KERNEL))
    stop ("'KERNEL' is missing in 'PredLik'")
  KERNEL <- match.fun(KERNEL)
  if(!any(names(formals(KERNEL))=="data"))
    stop ("'KERNEL' must have the argument 'data' in 'PredLik'")  
  if(!any(names(formals(KERNEL))=="log"))
    stop ("'KERNEL' must have the argument 'log' in 'PredLik'")  
   
  # marginal likelihood and scale from full sample
  ML.fs <- MargLik(N,mit=mit.fs,KERNEL,data=data.fs,...)
  scale.fs <-ML.fs$scale
  # marginal likelihood and scale from subsample
  ML.ss <- MargLik(N,mit=mit.ss,KERNEL,data=data.ss,...)
  scale.ss <-ML.ss$scale  
   
  # correct scaling
  scale.p <- scale.fs-scale.ss  
  # predictive likelihood from marginal likelihoods
  PL <- ML.fs$ML.mean/ML.ss$ML.mean
  
  return(list(PL=PL,scale=scale.p))
}
