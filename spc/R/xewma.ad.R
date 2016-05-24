# Computation of EWMA steady-state ARLs (mean monitoring)
xewma.ad <- function(l, c, mu1, mu0=0, zr=0, z0=0, sided="one", limits="fix", steady.state.mode="conditional", r=40) { 
  if ( l<=0 || l>1 ) 		stop("l has to be between 0 and 1")
  
  if ( c<=0 )			warning("usually, c has to be positive")
  
  if ( zr>c & sided=="one" )    stop("wrong reflexion border")
  
  if ( r<4 )			stop("r is too small")
  
  ctyp <- pmatch(sided, c("one", "two")) - 1
  if (is.na(ctyp))		stop("invalid ewma type")
  
  ltyp <- pmatch(limits, c("fix","vacl","fir","both","Steiner","stat")) - 1
  if ( is.na(ltyp) )		stop("invalid limits type")
  
  if ( (sided=="one") & !(limits %in% c("fix", "vacl", "stat")) )
				stop("not supported for one-sided EWMA (not reasonable or not implemented yet")

  styp <- pmatch(steady.state.mode, c("conditional", "cyclical")) - 1
  if (is.na(styp))		stop("invalid steady.state.mode")
  
  if ( abs(z0) > abs(c) ) 	stop("wrong restarting value")
  
  ad <- .C("xewma_ad", as.integer(ctyp), as.double(l),
           as.double(c), as.double(zr), as.double(mu0), as.double(mu1), as.double(z0),
           as.integer(ltyp), as.integer(styp), as.integer(r),
           ans=double(length=1), PACKAGE="spc")$ans 
  names(ad) <- "ad"
  return (ad)
}
