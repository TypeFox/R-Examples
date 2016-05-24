# Computation of CUSUM decision intervals -- alarm thresholds --  (variance monitoring)
scusum.crit <- function(k, L0, sigma, df, hs=0, sided="upper", mode="eq.tails", k2=NULL, hs2=0, r=40, qm=30) {
  if ( k<0 ) 
    stop("k has to be non-negative")
  if ( L0<1 )
    stop("L0 is too small")
  if ( hs<0 ) 
    stop("wrong headstart")
  if ( sided=="two" ) {
    if ( is.null(k2) )
      stop("in case of a two-sided CUSUM scheme one has to define two reference values")
    if ( k2<0 ) 
      stop("k2 has to be non-negative")
    if ( hs2<0 ) 
      stop("wrong headstart")
  }
    
  if ( sigma<=0 )
    stop("sigma must be positive")
  if ( df<1 )
    stop("df must be larger than or equal to 1")

  ctyp <- pmatch(sided, c("upper", "lower", "two")) - 1
  if ( is.na(ctyp) )
    stop("invalid cusum type")
  ltyp <- pmatch(mode, c("eq.tails", "unbiased")) - 1
  if ( is.na(ltyp) )
    stop("invalid limits cusum type")
  if ( r<10 ) 
    stop("r is too small")
  if ( qm<10 ) 
    stop("qm is too small")
   
  a.length <- 1
  if ( sided=="two" ) a.length <- 2
    
  h <- .C("scusum_crit", as.integer(ctyp),
          as.double(k), as.double(L0), as.double(hs),
          as.double(sigma), as.integer(df), as.integer(ltyp),
          as.double(k2), as.double(hs2),
          as.integer(r), as.integer(qm),
          ans=double(length=a.length), PACKAGE="spc")$ans
          
  if ( sided=="two" )  {
    names(h) <- c("hl","hu")
  } else {
    names(h) <- "h"
  }
  
  h
}
