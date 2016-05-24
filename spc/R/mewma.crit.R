# Computation of MEWMA threshold (multivariate mean monitoring)
mewma.crit <- function(l, L0, p, hs=0, r=20) {
  if ( l<=0 | l>1 )		stop("l has to be between 0 and 1")
  if ( L0<1 )			stop("L0 is too small")
  if ( p<1 )			stop("wrong dimension parameter")
  if ( hs<0 )			stop("wrong head start value")
  if ( r<4 )			stop("resolution too small")

  h <- .C("mewma_crit", as.double(l), as.double(L0),
                        as.integer(p), as.double(hs),
                        as.integer(r),
                        ans=double(length=1), PACKAGE="spc")$ans

  names(h) <- NULL
  h
}