# Computation of CUSUM ARLs (drift monitoring)
xDcusum.arl <- function(k, h, delta, hs=0, sided="one", mode="Gan", m=NULL, q=1, r=30, with0=FALSE) {
  if (k<0)
    stop("k has to be non-negative")
  if (h<=0)
    stop("h has to be positive")
  if ( hs<0 | (sided=="two" & hs>h/2+k) | (sided=="one" & hs>h/2+k) )
    stop("wrong headstart")
  if (r<4)
    stop("r is too small")
  if ( is.null(m) ) {
    m <- 0
  } else {
    if ( m<1 ) stop("m is too small") 
  }
  ctyp <- pmatch(sided, c("one", "two")) - 1
  if (is.na(ctyp))
    stop("invalid cusum type")
  cmode <- pmatch(mode, c("Gan", "Knoth")) - 1
  if (is.na(cmode))
    stop("invalid algorithm mode")
  q <- round(q)
  if (q<1)
    stop("wrong change point position (q)")
  arl <- .C("xDcusum_arl",as.integer(ctyp),as.double(k),
            as.double(h),as.double(hs),as.double(delta),as.integer(m),
            as.integer(r),as.integer(with0),as.integer(cmode),as.integer(q),
            ans=double(length=1),PACKAGE="spc")$ans
  names(arl) <- "arl"
  return (arl)
}
