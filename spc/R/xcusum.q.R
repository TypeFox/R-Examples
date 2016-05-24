# Computation of CUSUM quantiles (mean monitoring)
xcusum.q <- function(k, h, mu, alpha, hs=0, sided="one", r=40) {
  if ( k < 0 )  stop("k has to be non-negative")
  if ( h <= 0 ) stop("h has to be positive")
  if ( hs<0 | (sided=="two" & hs>h/2+k) | (sided=="one" & hs>h/2+k) ) stop("wrong headstart")
  if ( alpha <= 0 | alpha >= 1) stop("quantile level alpha must be in (0,1)")
  ctyp <- pmatch(sided, c("one", "two")) - 1
  if ( is.na(ctyp) ) stop("invalid cusum type")
  if ( r < 4 ) stop("r (dimension of Markov chain) is too small")
  quant <- .C("xcusum_q",
              as.integer(ctyp), as.double(k),as.double(h), as.double(alpha), as.double(hs), as.double(mu), as.integer(r),
              ans=double(length=1),PACKAGE="spc")$ans
  names(quant) <- "q"
  quant
}
