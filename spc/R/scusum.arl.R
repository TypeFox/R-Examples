# Computation of CUSUM ARLs (variance monitoring)
scusum.arl <- function(k, h, sigma, df, hs=0, sided="upper", k2=NULL, h2=NULL, hs2=0, r=40, qm=30, version=2) {
  if ( k<0 ) 
    stop("k has to be non-negative")
  if ( h<=0 ) 
    stop("h has to be positive")
  if ( hs<0 | hs>h ) 
    stop("wrong headstart")
  if ( sided=="two" ) {
    if ( is.null(k2) | is.null(h2) )
      stop("in case of a two-sided CUSUM scheme one has to define two sets of (k,h,hs)")
    if ( k2<0 ) 
      stop("k2 has to be non-negative")
    if ( h2<=0 ) 
      stop("h2 has to be positive")
    if ( hs2<0 | hs2>h2 ) 
      stop("wrong headstart")
  }
    
  if ( sigma<=0 )
    stop("sigma must be positive")
  if ( df<1 )
    stop("df must be larger than or equal to 1")

  ctyp <- pmatch(sided, c("upper", "lower", "two")) - 1
  if ( is.na(ctyp) )
    stop("invalid cusum type")
  if ( r<10 ) 
    stop("r is too small")
  if ( qm<10 ) 
    stop("qm is too small")

  arl <- .C("scusum_arl", as.integer(ctyp),
            as.double(k), as.double(h), as.double(hs),
            as.double(sigma), as.integer(df),
            as.double(k2), as.double(h2), as.double(hs2),
            as.integer(r), as.integer(qm), as.integer(version),
            ans=double(length=1), PACKAGE="spc")$ans 
  names(arl) <- "arl"
  
  arl
}
