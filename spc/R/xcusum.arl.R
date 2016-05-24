# Computation of CUSUM ARLs (mean monitoring)
xcusum.arl <- function(k, h, mu, hs=0, sided="one", method="igl", q=1, r=30) {
  if (k<0) 
    stop("k has to be non-negative")
  if (h<=0) 
    stop("h has to be positive")
  if ( hs<0 | (sided=="two" & hs>h/2+k) | (sided=="one" & hs>h) ) 
    stop("wrong headstart")
  if (r<4) 
    stop("r is too small")
  ctyp <- pmatch(sided, c("one", "two", "Crosier")) - 1
  if (is.na(ctyp)) 
    stop("invalid cusum type")
  mtyp <- pmatch(method, c("igl", "mc")) - 1
  if (is.na(mtyp)) 
    stop("invalid method")
  q <- round(q)
  if (q<1)
    stop("wrong change point position (q)")
  arl <- .C("xcusum_arl",
            as.integer(ctyp), as.double(k), as.double(h), as.double(hs), as.double(mu),
            as.integer(q), as.integer(r), as.integer(mtyp),
            ans=double(length=q), PACKAGE="spc")$ans
  names(arl) <- NULL
  return (arl)
}
