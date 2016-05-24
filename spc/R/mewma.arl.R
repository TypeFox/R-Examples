# Computation of MEWMA ARLs (multivariate mean monitoring)
mewma.arl <- function(l, cE, p, delta=0, hs=0, r=20, ntype=NULL, qm0=20, qm1=qm0) {
  if ( l<=0 | l>1 )		stop("l has to be between 0 and 1")
  if ( cE<=0 )			stop("threshold c has to be positive")
  if ( p<1 )			stop("wrong dimension parameter")
  if ( delta<0 )		stop("wrong magnitude value")
  if ( hs<0 )			stop("wrong head start value")
  if ( r<4 )			stop("resolution too small")
  if ( qm0<5 )			stop("more quadrature nodes needed")
  if ( qm1<5 )			stop("more quadrature nodes needed")
  
  if ( is.null(ntype) ) {
    if ( delta <1e-10 ) {
      ntype <- "gl2"
    } else {
      #if ( p %in% c(2,4) ) {
      if ( p==2 ) {
        ntype <- "gl3"
      } else {
        ntype <- "gl5"
      }
    }
  }
  
  qtyp <- pmatch(tolower(ntype), c("gl", "co", "ra", "cc", "mc", "sr", "co2", "gl2", "gl3", "gl4", "gl5", "co3", "co4")) - 1
  if ( is.na(qtyp) )		stop("invalid type of numerical algorithm")

  arl <- .C("mewma_arl", as.double(l), as.double(cE),
                         as.integer(p), as.double(delta),
                         as.double(hs), as.integer(r),
                         as.integer(qtyp), as.integer(qm0), as.integer(qm1),
                         ans=double(length=1), PACKAGE="spc")$ans

  names(arl) <- NULL
  arl
}