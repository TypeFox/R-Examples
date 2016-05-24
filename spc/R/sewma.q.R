# Computation of EWMA quantiles (variance monitoring)
sewma.q <- function(l, cl, cu, sigma, df, alpha, hs=1, sided="upper", r=40, qm=30) {
  if ( l <= 0 | l > 1 )			stop("l (lambda) has to be between 0 and 1")
  if ( cu<=0 )				stop("cu has to be positive")
  if ( cl<0 )				stop("cl has to be non-negative")
  if ( sided!="upper" & cl<1e-6 )	stop("cl is too small")  
  if ( sigma<=0 )			stop("sigma must be positive")
  if ( df<1 )				stop("df must be larger than or equal to 1")  
  if ( alpha <= 0 | alpha >= 1)		stop("quantile level alpha must be in (0,1)")
  if ( hs<cl | hs>cu )			stop("wrong headstart hs")  
  ctyp <- pmatch(sided, c("upper","Rupper","two","Rlower")) - 1
  if (is.na(ctyp))			stop("invalid ewma type") 
  if ( r<10 )				stop("r is too small")
  if ( qm<5 )				stop("qm is too small") 
  quant <- .C("sewma_q",
              as.integer(ctyp), as.double(l), as.double(cl), as.double(cu),
              as.double(alpha), as.double(hs), as.integer(r), as.double(sigma),
              as.integer(df), as.integer(qm),
              ans=double(length=1),PACKAGE="spc")$ans
  names(quant) <- "q"
  quant
}
