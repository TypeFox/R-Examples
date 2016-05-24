# Computation of EWMA critical values for given ARL 
# (simultaneous mean and variance monitoring)
xsewma.crit <- function(lx, ls, L0, df, mu0=0, sigma0=1, cu=NULL, hsx=0, hss=1,
                        s2.on=TRUE, sided="upper", mode="fixed",
                        Nx=30, Ns=40, qm=30) 
{
  if (lx<=0 || lx>1) 
    stop("lx has to be between 0 and 1")
  if (ls<=0 || ls>1) 
    stop("ls has to be between 0 and 1")
  if (L0<1) 
    stop("L0 is too small")
  if (sigma0<=0)
    stop("sigma0 must be positive")
  if (mode=="fixed" & sided=="two") {
    if (is.null(cu)) stop("set cu")
    if (cu<sigma0) stop("cu is too small")
    if (cu<=0) stop("cu must be positive")
    if (hss>cu) stop("hs must be smaller than cu")
    cu0 <- cu
  } else {
    cu0 <- 0
  }
  if (df<1)
    stop("df must be positive")
  s_squared <- as.numeric(s2.on)
  if ( !(s_squared %in% c(0,1)) )
    stop("wrong value for s2.on")
  ctyp <- pmatch(sided, c("upper","Rupper","two","lower")) - 1
  if (is.na(ctyp)) 
    stop("invalid ewma type")
  ltyp <- pmatch(mode, c("fixed","unbiased")) - 1
  if (is.na(ltyp)) 
    stop("invalid limits type")
  if (Nx<5) 
    stop("r.x is too small")
  if (Ns<10) 
    stop("r.s is too small")
  if (qm<10) 
    stop("qm is too small")
  c <- .C("xsewma_crit",as.integer(ctyp),as.integer(ltyp),
          as.double(lx),as.double(ls),
          as.double(L0),as.double(cu0),as.double(hsx),as.double(hss),
          as.double(mu0),as.double(sigma0),
          as.integer(df),as.integer(Nx),as.integer(Ns),
          as.integer(qm),
          ans=double(length=3),PACKAGE="spc")$ans
  names(c) <- c("cx","cl","cu")
  return (c)
}

