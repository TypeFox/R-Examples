# Computation of EWMA survival function (mean monitoring)
xewma.sf <- function(l, c, mu, n, zr=0, hs=0, sided="one", limits="fix", q=1, r=40) {
  if ( l <= 0 | l > 1 )          stop("l (lambda) has to be between 0 and 1")
  if ( c <= 0 )                  warning("usually, c has to be positive")
  if ( n < 1 )                   stop("n has to be a natural number")
  if ( zr > c & sided == "one")  stop("wrong reflexion border")
  if ( (sided == "two" & abs(hs) > c) | (sided == "one" & ( hs < zr | hs > c )) )
                                 warning("unusual headstart")
  ctyp <- pmatch(sided, c("one", "two")) - 1
  if ( is.na(ctyp) )             stop("invalid ewma type")
  ltyp <- -1 + pmatch(limits, c("fix", "vacl", "fir", "both", "Steiner", "stat", "test"))
  if (is.na(ltyp))               stop("invalid limits type")
  if ( (sided=="one") & !( limits %in% c("fix", "vacl", "stat") ) )
                                 stop("not supported for one-sided EWMA (not reasonable or not implemented yet")
  if ( r < 4 )                   stop("r is too small")
  q <- round(q)
  if ( q<1 )                     stop("wrong change point position (q)")
  sf <- .C("xewma_sf",
           as.integer(ctyp), as.double(l), as.double(c), as.double(zr), as.double(hs), as.double(mu),
           as.integer(ltyp), as.integer(r), as.integer(q), as.integer(n),
           ans=double(length=n),PACKAGE="spc")$ans
  names(sf) <- NULL
  sf
}
