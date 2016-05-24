# Computation of res-EWMA ARLs (scale monitoring)
s.res.ewma.arl <- function(l,cu,sigma,mu=0,alpha=0,n=5,hs=1,r=40,qm=30) {
  if ( l <= 0 || l > 1 )
    stop("l has to be between 0 and 1")
  if ( cu <= 0 )
    warning("usually, cu has to be positive")
  if ( sigma <= 0 )
    stop("sigma must be positive")
  if ( abs(alpha) > 1 )
    warning("nonstationary AR(1) process")
  if ( n < 2 )
    warning("n is too small")
  n <- round(n)
  if ( abs(hs) > cu ) 
    warning("unusual headstart")
  if ( r < 4 )
    stop("r is too small")
  if ( qm < 10 ) 
    stop("qm is too small")
  ctyp <- 1 # later more
  arl <- .C("s_res_ewma_arl",as.double(alpha),as.integer(n-1),
            as.integer(ctyp),as.double(l),
            as.double(cu),as.double(hs),
            as.double(sigma),as.double(mu),as.integer(r),as.integer(qm),
            ans=double(length=1),PACKAGE="spc")$ans 
  names(arl) <- "arl"
  return (arl)
}

