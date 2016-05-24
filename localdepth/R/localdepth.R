#############################################################
#
#	localdepth function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: November, 06, 2013
#	Version: 0.4
#
#	Copyright (C) 2013 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

localdepth <- function(x, y=NULL, tau, use=c('volume', 'diameter'), method=c('simplicial', 'ellipsoid', 'halfspace', 'mahalanobis', "hyperspheresimplicial"), type=c('exact', 'approx'), nsamp='all', nmax=1, tol=10^(-9), dimension=NULL, location=NULL, covariance=NULL) {
  use <- match.arg(use)
  method <- match.arg(method)
  type <- match.arg(type)

## simplicial
  if (method=='simplicial') {
    if (is.circular(x))
      localdepth.simp.circular(x=x, y=y, tau=tau, use=use)
    else if (type=='exact') {
      if (NCOL(x) < 3 & nsamp=='all') 
        localdepth.simp(x=x, y=y, tau=tau, use=use)
      else
        localdepth.simp.exact(x=x, y=y, tau=tau, use=use, nsamp=nsamp, nmax=nmax, tol=tol)
    } else
      localdepth.simp.approx(x=x, y=y, tau=tau, use=use, nsamp=nsamp, nmax=nmax, tol=tol)
## halfspace    
  } else if (method=='halfspace') {
      if (is.circular(x))
        localdepth.tukey.circular(x=x, y=y, tau=tau, tol=tol)
      else 
        localdepth.halfspace(x=x, y=y, tau=tau) 
## mahalanobis    
  } else if (method=='mahalanobis') {
      if (is.circular(x)) stop("method 'mahalanobis' is not implemented for circular data")
      localdepth.mahalanobis(x=x, y=y, tau=tau, nsamp=nsamp, nmax=nmax, location=location, covariance=covariance)
## ellipsoid
  } else if (method=='ellipsoid')  {
      if (is.circular(x)) stop("method 'ellipsoid' is not implemented for circular data")
      localdepth.ellipsoid(x=x, y=y, tau=tau, use=use, nsamp=nsamp, nmax=nmax, tol=tol, dimension=dimension)
  } else  if (method=='halfspace')  {
    if (is.circular(x))
      localdepth.tukey.circular(x, y=y, tau=tau, tol=tol)
    if (type=='exact')
      localdepth.halfspace(x=x, y=y, tau=tau, use=use)
    else
      stop("method 'halfspace' is not implemented for approximated calculation")
  } else {
    if (is.circular(x))
      localdepth.simp.circular(x=x, y=y, tau=tau, use=use)
    if (is.null(dim(x))) {
      warning("method='hyperspheresimplicial' is available only for vector on dimension 2 or larger, we use 'simplicial' instead and nsamp is set to 'all'")
      localdepth.simp(x=x, y=y, tau=tau, use=use)
    } else if (type=='exact')
      localdepth.simp.exact.hyperspheres(x=x, y=y, tau=tau, use=use, nsamp=nsamp, nmax=nmax, tol=tol)
    else
      stop("method 'hyperspheresimplicial' is not implemented for approximated calculation")

  }
}

