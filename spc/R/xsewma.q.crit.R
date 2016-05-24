# Computation of EWMA critical values for given QRL 
# (simultaneous mean and variance monitoring)
xsewma.q.crit <- function(lx, ls, L0, alpha, df, mu0=0, sigma0=1, csu=NULL, hsx=0, hss=1, sided="upper", mode="fixed", Nx=20, Ns=40, qm=30, c.error=1e-12, a.error=1e-9)
{
  if (lx<=0 || lx>1)			stop("lx has to be between 0 and 1")
  if (ls<=0 || ls>1)			stop("ls has to be between 0 and 1")
  if (L0<1)				stop("L0 is too small")
  if ( alpha<=0 | alpha>=1 )		stop("quantile level alpha must be in (0,1)")
  if ( df<1 )				stop("df must be positive")
  if ( sigma0<=0 )			stop("sigma0 must be positive")  
  if ( mode=="fixed" & sided=="two" ) {
    if ( is.null(csu) )			stop("set csu")
    if ( csu<sigma0 )			stop("csu is too small")
    if ( csu<=0 )			stop("csu must be positive")
    if ( hss>csu )			stop("hs must be smaller than csu")
    cu0 <- csu
  } else {
    cu0 <- 0
  }
  ctyp <- pmatch(sided, c("upper","two")) - 1
  if ( is.na(ctyp) )			stop("invalid ewma type")
  ltyp <- pmatch(mode, c("fixed","unbiased")) - 1
  if ( is.na(ltyp) )			stop("invalid limits type")
  if ( Nx<5 )				stop("Nx is too small")
  if ( Ns<10 )				stop("Ns is too small")
  if ( qm<10 )				stop("qm is too small")
  c <- .C("xsewma_q_crit", as.integer(ctyp), as.integer(ltyp),
          as.double(lx), as.double(ls),
          as.double(L0), as.double(alpha),
          as.double(cu0), as.double(hsx), as.double(hss),
          as.double(mu0), as.double(sigma0), as.integer(df), as.integer(Nx), as.integer(Ns), as.integer(qm),
          as.double(c.error), as.double(a.error),
          ans=double(length=3),PACKAGE="spc")$ans
  names(c) <- c("cx", "csl","csu")
  return (c)
}

