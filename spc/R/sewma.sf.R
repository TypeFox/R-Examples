# Computation of EWMA survival function (variance monitoring)
sewma.sf <- function(n, l, cl, cu, sigma, df, hs=1, sided="upper", r=40, qm=30) {
  if ( n < 1 )			stop("n has to be a natural number")
  if ( l <= 0 | l > 1 )		stop("l (lambda) has to be between 0 and 1")
  if ( cu<=0 )			stop("cu has to be positive")
  if ( cl<0 )			stop("cl has to be non-negative")
  if ( sided!="upper" & cl<1e-6 ) stop("cl is too small")
  if ( sigma<=0 )		stop("sigma must be positive")
  if ( df<1 )			stop("df must be larger than or equal to 1")  
  if ( hs<cl | hs>cu )		stop("wrong headstart hs")
  if ( r<10 )			stop("r is too small")
  ctyp <- pmatch(sided, c("upper","Rupper","two","Rlower")) - 1
  if (is.na(ctyp))		stop("invalid ewma type") 
  if ( qm<5 )			stop("qm is too small") 
  sf <- .C("sewma_sf",
           as.integer(ctyp), as.double(l), as.double(cl), as.double(cu), as.double(hs), as.integer(r),
           as.double(sigma), as.integer(df), as.integer(qm), as.integer(n),
           ans=double(length=n),PACKAGE="spc")$ans
  names(sf) <- NULL
  sf
}
