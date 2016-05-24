#############################################################
#	All functions in this file are copyrighted
#	Author: C. Agostinelli, A. Marazzi,
#               V.J. Yohai and A. Randriamiharisoa
#	Maintainer e-mail: claudio@unive.it
#	Date: January, 1, 2013
#	Version: 0.1
#	Copyright (C) 2013 C. Agostinelli A. Marazzi,
#                  V.J. Yohai and A. Randriamiharisoa
#############################################################

"Regtau.f" <- function(x,y,b1,c1,b2,c2,N,tol=1e-6,seed=567) {
  n  <-length(y)
  rs <- tmp1 <- tmp2 <- rep(0,n)
  sa <- sb <- ta <- rep(0,N)
  f.res <- .Fortran("regtau",
    x=as.double(x),
    y=as.double(y),
    n=to.integer(n),
    nn=to.integer(N),
    b1=to.single(b1),
    c1=to.single(c1),
    b2=to.single(b2),
    c2=to.single(c2),
    tol=to.single(tol),
    iseed=to.integer(seed),
    ao=double(1),
    bo=double(1),
    to=double(1),
    rs=as.double(rs),
    sa=as.double(sa),
    sb=as.double(sb),
    ta=as.double(ta),
    tmp1=to.single(tmp1),
    tmp2=to.single(tmp2),
    PACKAGE = "robustloggamma")
  result <- list(ao=f.res$ao,bo=f.res$bo,to=f.res$to)
  return(result)  
}


