# Computation of EWMA critical values for given QRL (variance monitoring)
sewma.q.crit <- function(l, L0, alpha, df, sigma0=1, cl=NULL, cu=NULL, hs=1, sided="upper", mode="fixed", ur=4, r=40, qm=30, c.error=1e-12, a.error=1e-9) {
  cu0 <- cl0 <- 0
  if ( l<=0 | l>1 )		stop("l has to be between 0 and 1")
  if ( L0<1 )			stop("L0 is too small")
  if ( alpha<=0 | alpha>=1 )	stop("quantile level alpha must be in (0,1)")
  if ( df<1 )			stop("df must be positive")
  if ( sigma0<=0 )		stop("sigma0 must be positive")
  if ( sided=="Rupper" ) {
    if ( is.null(cl) )		stop("set cl")
    if ( cl<=0 )		stop("cl must be positive")
    cl0 <- cl
    if ( hs<cl0 )		stop("hs must not be smaller than cl")
  }
  if ( sided=="Rlower" ) {
    if ( is.null(cu) )		stop("set cu")
    if ( cu<sigma0 )		stop(paste("cu must be larger than sigma0 =", sigma0))
    cu0 <- cu
    if ( hs>cu0 )		stop("hs must not be larger than cu")
  }
  if ( sided=="two" & mode=="fixed" ) {
    if ( is.null(cu) )		stop("set cu")
    if ( cu<sigma0 )		stop(paste("cu must be larger than sigma0 =", sigma0))
    cu0 <- cu
    if ( hs>cu0 )		stop("hs must not be larger than cu")
  }
  ctyp <- pmatch(sided, c("upper", "Rupper", "two", "Rlower")) - 1
  if ( is.na(ctyp) )		stop("invalid ewma type")
  ltyp <- pmatch(mode, c("fixed", "unbiased", "classic")) - 1
  if ( is.na(ltyp) )		stop("invalid limits type")
  if ( r<10 )			stop("r is too small")
  if ( qm<10 )			stop("qm is too small")
  c <- .C("sewma_q_crit",
          as.integer(ctyp), as.integer(ltyp), as.double(l),
          as.integer(L0), as.double(alpha),
          as.double(cl0), as.double(cu0), as.double(hs),
          as.double(sigma0), as.integer(df), as.integer(r), as.integer(qm), as.double(ur),
          as.double(c.error), as.double(a.error),
          ans=double(length=2),PACKAGE="spc")$ans
  names(c) <- c("cl", "cu")
  return (c)
}

