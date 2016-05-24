# Computation of EWMA ARLs (variance monitoring) with pre-run uncertainty
sewma.arl.prerun <- function(l, cl, cu, sigma, df1, df2, hs=1, sided="upper", r=40, qm=30, qm.sigma=30, truncate=1e-10) {
  if ( l<=0 || l>1 ) 
    stop("l has to be between 0 and 1")
  if ( cl<0 )
    stop("cl has to be non-negative")
  if ( cu<=0 ) 
    stop("cu has to be positive")
  if ( sigma<=0 )
    stop("sigma must be positive")
  if ( df1<1 )
    stop("df1 must be larger than or equal to 1")
  if ( df2<1 )
    stop("df2 must be larger than or equal to 1")
  if ( hs<cl | hs>cu ) 
    stop("wrong headstart")
  ctyp <- pmatch(sided, c("upper", "Rupper", "two", "Rlower")) - 1
  if (is.na(ctyp))
    stop("invalid ewma type")
  if ( r<10 ) 
    stop("r is too small")
  if ( qm<10 ) 
    stop("qm is too small")
  if ( qm.sigma<4 )
    stop("qm.sigma is too small")
  if ( truncate < 0 | truncate >= 0.5 )
    stop("wrong value for truncate (should follow 0 < truncate < 0.5)")  
  arl <- .C("sewma_arl_prerun", as.integer(ctyp), as.double(l),
            as.double(cl), as.double(cu), as.double(hs),
            as.double(sigma), as.integer(df1), as.integer(r), as.integer(qm), 
            as.integer(df2), as.integer(qm.sigma), as.double(truncate),
            ans=double(length=1),PACKAGE="spc")$ans 
  names(arl) <- "arl"
  return (arl)
}
