#############################################################
#
#	localdepth.similarity function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: April, 02, 2009
#	Version: 0.2-1
#
#	Copyright (C) 2009 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

localdepth.similarity <- function(x, y=NULL, tau, use=c('volume', 'diameter'), method=c('simplicial', 'ellipsoid', 'mahalanobis'), type=c('exact', 'approx'), nsamp='all', nmax=1, tol=10^(-9), dimension=NULL, location=NULL, covariance=NULL, weight=NULL) {
  use <- match.arg(use)
  method <- match.arg(method)
  type <- match.arg(type)

## simplicial
  if (method=='simplicial') {
    if (is.circular(x))
      localdepth.similarity.simp.circular(x=x, y=y, tau=tau, use=use)
    if (type=='exact') {
      if (NCOL(x) < 3 & nsamp=='all')
        localdepth.similarity.simp(x=x, y=y, tau=tau, use=use, weight=weight)
      else
        localdepth.similarity.simp.exact(x=x, y=y, tau=tau, use=use, nsamp=nsamp, nmax=nmax, tol=tol)
    } else
      localdepth.similarity.simp.approx(x=x, y=y, tau=tau, use=use, nsamp=nsamp, nmax=nmax, tol=tol)
## mahalanobis    
  } else if (method=='mahalanobis') {
      if (is.circular(x)) stop("method 'mahalanobis' is not implemented for circular data")
      localdepth.similarity.mahalanobis(x=x, y=y, tau=tau, nsamp=nsamp, nmax=nmax, location=location, covariance=covariance, weight=weight)
## ellipsoid
  } else {
      if (is.circular(x)) stop("method 'ellipsoid' is not implemented for circular data")
      localdepth.similarity.ellipsoid(x=x, y=y, tau=tau, use=use, nsamp=nsamp, nmax=nmax, tol=tol, dimension=dimension)
  }
}

