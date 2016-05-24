# Computation of EWMA critical values for given ARL (variance monitoring) based on ln S^2
lns2ewma.crit <- function(l, L0, df, sigma0=1, cl=NULL, cu=NULL, hs=NULL, sided="upper", mode="fixed", r=40) {

  #mitte <- -1/df - 1/3/df^2 + 2/15/df^4 # approx following Crowder/Hamilton
  mitte <- log(2/df) + digamma(df/2)
  
  if ( is.null(hs) ) hs <- mitte
  cu0 <- cl0 <- 0
  
  if ( l<=0 || l>1 ) stop("l has to be between 0 and 1")
  
  if ( L0<1 ) stop("L0 is too small")
  
  if ( df<1 ) stop("df must be positive")
  
  if ( sigma0<=0 ) stop("sigma0 must be positive")
    
  if ( sided=="upper" ) {
    if ( is.null(cl) ) cl <- mitte
    #if ( cl > mitte + 1e-9 ) stop(paste("cl has to be smaller than", mitte))
    cl0 <- cl
    if ( hs<cl0-1e-9 ) stop("hs must not be smaller than cl")
  }
  if ( sided=="lower" ) {
    if ( is.null(cu) ) cu <- mitte
    #if ( cu < mitte - 1e-9 ) stop(paste("cu has to be larger than", mitte))
    cu0 <- cu
    if ( hs>cu0+1e-9 ) stop("hs must not be larger than cu")
  }
  if (sided=="two" & mode=="fixed") {
    if ( is.null(cu) ) stop("set cu")
    #if ( cu < mitte - 1e-9 ) stop(paste("cu has to be larger than", mitte))
    cu0 <- cu
    if ( hs>cu0+1e-9 ) stop("hs must not be larger than cu")
  }
  
  ctyp <- pmatch(sided, c("upper", "lower", "two")) - 1
  if (is.na(ctyp)) stop("invalid ewma type")
  
  ltyp <- pmatch(mode, c("fixed", "unbiased", "eq.tails", "vanilla")) - 1
  if (is.na(ltyp)) stop("invalid limits type")
  
  if ( r<10 ) stop("r is too small")
  
  c <- .C("lns2ewma_crit", as.integer(ctyp), as.integer(ltyp), as.double(l),
          as.double(L0), as.double(cl0), as.double(cu0), as.double(hs),
          as.double(sigma0), as.integer(df), as.integer(r),
          ans=double(length=2),PACKAGE="spc")$ans
          
  names(c) <- c("cl", "cu")
  return (c)
}

