# Computation of res-EWMA ARLs (mean monitoring)
x.res.ewma.arl <- function(l, c, mu, alpha=0, n=5, hs=0, r=40) {
  if ( l <= 0 || l > 1 )
    stop("l has to be between 0 and 1")
  if ( c <= 0 )
    warning("usually, c has to be positive")
  if ( abs(alpha) > 1 )
    warning("nonstationary AR(1) process")
  if ( n < 1 )
    warning("n is too small")
  n <- round(n)
  if ( abs(hs) > c ) 
    warning("unusual headstart")
  if ( r < 4 )
    stop("r is too small")
  ctyp <- 1 # later more
  arl <- .C("x_res_ewma_arl",as.double(alpha),as.integer(n),
            as.integer(ctyp),as.double(l),
            as.double(c),as.double(hs),
            as.double(mu),as.integer(r),
            ans=double(length=1),PACKAGE="spc")$ans 
  names(arl) <- "arl"
  return (arl)
}