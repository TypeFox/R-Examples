# Computation of EWMA ARLs (simultaneous mean & variance monitoring)
xsewma.arl <- function(lx, cx, ls, csu, df, mu, sigma,
                       hsx=0, Nx=40,
                       csl=0, hss=1, Ns=40,
                       s2.on=TRUE, sided="upper", qm=30) {
  if (lx<=0 | lx>1) 
    stop("lx has to be between 0 and 1")
  if (ls<=0 | ls>1)
    stop("ls has to be between 0 and 1")
  if (cx<=0)
    stop("cx has to be positive")
  if (csu<=0) 
    stop("csu has to be positive")
  if (csl<0)
    stop("clu has to be non-negative")
  if ( sigma<=0 )
    stop("sigma must be positive")
  if ( df<1 )
    stop("df must be larger than or equal to 1")
  s_squared <- as.numeric(s2.on)
  if ( !(s_squared %in% c(0,1)) )
    stop("wrong value for s2.on")
  if ( abs(hsx)>cx )
    stop("wrong headstart hsx")
  if ( hss<csl | hss>csu ) 
    stop("wrong headstart hss")
  if (Nx<5)
    stop("Nx is too small")
  if (Ns<10) 
    stop("Ns is too small")
  if (qm<5) 
    stop("qm is too small")
  ctyp <- pmatch(sided, c("upper","Rupper","two","lower")) - 1
  if (is.na(ctyp)) 
    stop("invalid ewma type")
  arl <- .C("xsewma_arl",as.integer(ctyp),
            as.double(lx),as.double(cx),as.double(hsx),as.integer(Nx),
            as.double(ls),as.double(csl),as.double(csu),as.double(hss),
            as.integer(Ns),
            as.double(mu),as.double(sigma),
            as.integer(df),as.integer(qm),
            as.integer(s_squared),
            ans=double(length=1),PACKAGE="spc")$ans 
  names(arl) <- "arl"
  return (arl)
}
