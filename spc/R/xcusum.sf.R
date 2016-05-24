# Computation of CUSUM survival function (mean monitoring)
xcusum.sf <- function(k, h, mu, n, hs=0, sided="one", r=40) {
  if ( k < 0 )  stop("k has to be non-negative")
  if ( h <= 0 ) stop("h has to be positive")
  if ( hs < 0 | (sided=="two" & hs>h/2+k) | (sided=="one" & hs>h/2+k) ) stop("wrong headstart")
  if ( n < 1 ) stop("n has to be a natural number")
  ctyp <- pmatch(sided, c("one", "two")) - 1
  if ( is.na(ctyp) ) stop("invalid cusum type")
  if ( r < 4 ) stop("r is too small")
  sf <- .C("xcusum_sf",
   as.integer(ctyp), as.double(k), as.double(h), as.double(hs), as.double(mu), as.integer(r), as.integer(n),
              ans=double(length=n),PACKAGE="spc")$ans
  names(sf) <- NULL
  sf
}
