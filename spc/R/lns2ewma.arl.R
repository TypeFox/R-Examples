# Computation of EWMA ARLs (variance monitoring) based on ln S^2
lns2ewma.arl <- function(l, cl, cu, sigma, df, hs=NULL, sided="upper", r=40) {

  #mitte <- -1/df - 1/3/df^2 + 2/15/df^4 # approx following Crowder/Hamilton
  mitte <- log(2/df) + digamma(df/2)
  
  if ( is.null(cl) ) cl <- mitte
  
  if ( is.null(cu) ) cu <- mitte
  
  if ( is.null(hs) ) hs <- mitte
  
  if ( l<=0 || l>1 ) stop("l has to be between 0 and 1")
  
  #if ( cu < mitte )  stop(paste("cu has to be larger than", mitte))
  #if ( cl > mitte )  stop(paste("cl has to be smaller than", mitte))
  
  if ( sigma<=0 ) stop("sigma must be positive")
  
  if ( df<1 ) stop("df must be larger than or equal to 1")
  
  if ( hs<cl-1e-9 | hs>cu+1e-9 ) stop("wrong headstart")
  
  ctyp <- pmatch(sided, c("upper", "lower", "two")) - 1
  if ( is.na(ctyp) ) stop("invalid ewma type")
  
  if ( r<10 ) stop("r is too small")
  
  arl <- .C("lns2ewma_arl", as.integer(ctyp), as.double(l),
            as.double(cl), as.double(cu), as.double(hs),
            as.double(sigma), as.integer(df), as.integer(r),
            ans=double(length=1),PACKAGE="spc")$ans 
            
  names(arl) <- "arl"
  return (arl)
}
