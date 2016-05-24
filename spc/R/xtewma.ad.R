# Computation of EWMA steady-state ARLs (mean monitoring, t distributed data)
xtewma.ad <- function(l, c, df, mu1, mu0=0, zr=0, z0=0, sided="one", limits="fix", steady.state.mode="conditional", mode="tan", r=40) { 
  if ( l<=0 || l>1 ) 		warning("l has to be between 0 and 1")
  
  if ( c<=0 )			warning("usually, c has to be positive")
  
  if ( zr>c & sided=="one" )    stop("wrong reflexion border")
  
  if ( r<4 )			stop("r is too small")
  
  ctyp <- pmatch(sided, c("one", "two")) - 1
  if (is.na(ctyp))		stop("invalid ewma type")
  
  ltyp <- pmatch(limits, c("fix","vacl")) - 1
  if ( is.na(ltyp) )		stop("invalid limits type")


  styp <- pmatch(steady.state.mode, c("conditional", "cyclical")) - 1
  if (is.na(styp))		stop("invalid steady.state.mode")
  
  if ( abs(z0) > abs(c) ) 	stop("wrong restarting value")
  
  ntyp <- -1 + pmatch(mode, c("identity", "sin", "sinh", "tan"))
  if ( is.na(ntyp) )		stop("substitution type not provided (yet)")    
  
  ad <- .C("xtewma_ad", as.integer(ctyp), as.double(l),
           as.double(c), as.double(zr), as.integer(df), as.double(mu0), as.double(mu1), as.double(z0),
           as.integer(ltyp), as.integer(styp), as.integer(r), as.integer(ntyp),
           ans=double(length=1), PACKAGE="spc")$ans 
  names(ad) <- "ad"
  return (ad)
}
