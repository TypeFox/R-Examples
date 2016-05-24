# Computation of res-EWMA ARLs (simultaneous mean & variance monitoring)
xs.res.ewma.arl <- function(lx, cx, ls, csu, mu, sigma, alpha=0, n=5,
                       hsx=0, rx=40, hss=1, rs=40, qm=30) {
  if ( lx<=0 || lx>1 ) 
    stop("lx has to be between 0 and 1")
  if ( ls<=0 || ls>1 )
    stop("ls has to be between 0 and 1")
  if ( cx <= 0 )
    stop("cx has to be positive")
  if ( csu <= 0 ) 
    stop("csu has to be positive")
  if ( sigma <= 0 )
    stop("sigma must be positive")
  if ( abs(alpha)>1 )
    warning("nonstationary AR(1) process")
  if ( n < 2 )
    warning("n is too small")
  n <- round(n)
  if ( abs(hsx) > cx )
    stop("wrong headstart hsx")
  if ( hss < 0 | hss > csu ) 
    stop("wrong headstart hss")
  if ( rx < 5 )
    stop("rx is too small")
  if ( rs < 10 ) 
    stop("rs is too small")
  if ( qm < 5 ) 
    stop("qm is too small")
  ctyp <- 1 # later more
  arl <- .C("xsewma_res_arl",as.double(alpha),as.integer(n-1),as.integer(ctyp),
            as.double(lx),as.double(cx),as.double(hsx),as.integer(rx),
            as.double(ls),as.double(csu),as.double(hss),as.integer(rs),
            as.double(mu),as.double(sigma),as.integer(qm),
            ans=double(length=1),PACKAGE="spc")$ans 
  names(arl) <- "arl"
  return (arl)
}
