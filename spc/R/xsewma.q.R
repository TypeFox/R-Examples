# Computation of EWMA RL quantiles (simultaneous mean & variance monitoring)
xsewma.q <- function(lx, cx, ls, csu, df, alpha, mu, sigma, hsx=0, Nx=40, csl=0, hss=1, Ns=40, sided="upper", qm=30) {
  if ( lx<=0 | lx>1 ) 
    stop("lx has to be between 0 and 1")
  if ( ls<=0 | ls>1 )
    stop("ls has to be between 0 and 1")
  if ( cx<=0 )
    stop("cx has to be positive")
  if ( csu<=0 ) 
    stop("csu has to be positive")
  if ( df<1 )
    stop("df must be larger than or equal to 1")
  if ( alpha <= 0 | alpha >= 1)
    stop("quantile level alpha must be in (0,1)") 
  if ( sigma<=0 )
    stop("sigma must be positive")
  if ( abs(hsx)>cx )
    stop("wrong headstart hsx")
  if ( Nx<5 )
    stop("Nx is too small")    
  if ( csl<0 )
    stop("clu has to be non-negative")  
  if ( hss<csl | hss>csu ) 
    stop("wrong headstart hss")
  if ( Ns<10 ) 
    stop("Ns is too small")
  ctyp <- pmatch(sided, c("upper","two")) - 1
  if (is.na(ctyp)) 
    stop("invalid ewma type")
  if ( qm<5 ) 
    stop("qm is too small")
  quant <- .C("xsewma_q",as.integer(ctyp),as.double(alpha),
            as.double(lx),as.double(cx),as.double(hsx),as.integer(Nx),
            as.double(ls),as.double(csl),as.double(csu),as.double(hss),as.integer(Ns),
            as.double(mu),as.double(sigma),
            as.integer(df),as.integer(qm),
            ans=double(length=1),PACKAGE="spc")$ans 
  names(quant) <- "q"
  quant
}
