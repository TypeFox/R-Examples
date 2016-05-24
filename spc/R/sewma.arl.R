# Computation of EWMA ARLs (variance monitoring)
sewma.arl <- function(l, cl, cu, sigma, df, s2.on=TRUE, hs=NULL, sided="upper", r=40, qm=30) {
  mitte <- sqrt( 2/df ) * gamma( (df+1)/2 )/ gamma( df/2 )
  if ( is.null(hs) ) {
    if ( s2.on ) { hs <- 1 } else { hs <- mitte }
  }
  
  if ( l<=0 || l>1 ) 
    stop("l has to be between 0 and 1")
  if ( cu<=0 ) 
    stop("cu has to be positive")
  if ( cl<0 )
    stop("cl has to be non-negative")
  if ( sigma<=0 )
    stop("sigma must be positive")
  if ( df<1 )
    stop("df must be larger than or equal to 1")
  s_squared <- as.numeric(s2.on)
  if ( !(s_squared %in% c(0,1)) )
    stop("wrong value for s2.on")
  
  if ( hs<cl | hs>cu ) 
    stop("wrong headstart")
  ctyp <- pmatch(sided, c("upper", "Rupper", "two", "Rlower")) - 1
  if (is.na(ctyp))
    stop("invalid ewma type")
  if (r<10) 
    stop("r is too small")
  if (qm<10) 
    stop("qm is too small")
  arl <- .C("sewma_arl",as.integer(ctyp),as.double(l),
            as.double(cl),as.double(cu),as.double(hs),
            as.double(sigma),as.integer(df),as.integer(r),as.integer(qm),
            as.integer(s_squared),
            ans=double(length=1),PACKAGE="spc")$ans 
  names(arl) <- "arl"
  return (arl)
}
