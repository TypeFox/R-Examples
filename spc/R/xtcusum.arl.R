# Computation of CUSUM ARLs (mean monitoring, t distributed data)
xtcusum.arl <- function(k, h, df, mu, hs=0, sided="one", mode="tan", r=30) {
  if ( k < 0 ) 
    stop("k has to be non-negative")
  if ( h <= 0 ) 
    stop("h has to be positive")
  if ( df < 1 )
    stop("df must be greater or equal to 1")
  if ( hs < 0 | ( sided=="two" & hs>h/2+k ) | ( sided=="one" & hs>h ) ) 
    stop("wrong headstart")
    ntyp <- -1 + pmatch(mode, c("identity", "sin", "sinh", "tan"))
  if ( is.na(ntyp) )
    stop("substitution type not provided (yet)")    
  if ( r < 4 ) 
    stop("r is too small")
  ctyp <- pmatch(sided, c("one", "two")) - 1
  if (is.na(ctyp)) 
    stop("invalid cusum type")
  arl <- .C("xtcusum_arl",
            as.integer(ctyp), as.double(k), as.double(h), as.double(hs), as.integer(df), double(mu), as.integer(r), as.integer(ntyp),
            ans=double(length=1), PACKAGE="spc")$ans
  names(arl) <- NULL
  return (arl)
}