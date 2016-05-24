#############################################################
#
#	quantile.localdepth function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: November, 06, 2013
#	Version: 0.2
#
#	Copyright (C) 2013 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

quantile.localdepth <- function(x, probs, use=c('volume', 'diameter'),  method=c('simplicial', 'ellipsoid', 'halfspace', 'mahalanobis', 'hyperspheresimplicial'), nsamp='all', size=FALSE, dimension=NULL, covariance=NULL, ...) {
  use <- match.arg(use)
  method <- match.arg(method)
## simplicial and halfspace
  if (method=='simplicial' | method=='halfspace') {
    if (is.circular(x))
      quantile.simp.circular(x=x, probs=probs, all=size, ...)
    else
      quantile.simp(x=x, probs=probs, use=use, nsamp=nsamp, all=size, ...)
## mahalanobis    
  } else if (method=='mahalanobis') {
      if (is.circular(x)) stop("method 'mahalanobis' is not implemented for circular data")
      quantile.mahalanobis(x=x, probs=probs, nsamp=nsamp, all=size, covariance=covariance, ...)
## ellipsoid
  } else if (method=='mahalanobis') {
      if (is.circular(x)) stop("method 'ellipsoid' is not implemented for circular data")
      quantile.ellipsoid(x=x, probs=probs, use=use, nsamp=nsamp, all=size, dimension=dimension, ...)
  } else {
      if (is.circular(x))
        quantile.simp.circular(x=x, probs=probs, all=size, ...)        
      else
        quantile.hyperspheres(x=x, probs=probs, use=use, nsamp=nsamp, size=size, ...)
  }
}
