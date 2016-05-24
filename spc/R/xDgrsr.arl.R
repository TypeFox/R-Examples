# Computation of GRSR (Girshick, Rubin, Shiryaev, Roberts) ARLs (drift monitoring)
xDgrsr.arl <- function(k, g, delta, zr=0, hs=NULL, sided="one", m=NULL, mode="Gan", q=1, r=30, with0=FALSE) {
  if (k<0)
    stop("k has to be non-negative")
  if (g<=0)
    stop("g has to be positive")
  if (zr>g)
    stop("zr has to be smaller than g")
  if ( !is.null(hs) ) {
    if ( hs>g )
      stop("wrong headstart")
  } else {
    hs <- 2*g # mimics hs = -inf
  }
  ctyp <- pmatch(sided, c("one", "two")) - 1
  if (is.na(ctyp))
    stop("invalid grsr type")
  if (r<4)
    stop("r is too small")
  if ( is.null(m) ) {
    m <- 0
  } else {
    if ( m<1 ) stop("m is too small") 
  }
  cmode <- pmatch(mode, c("Gan", "Knoth")) - 1
  if (is.na(cmode))
    stop("invalid algorithm mode")
  q <- round(q)
  if (q<1)
    stop("wrong change point position (q)")
  arl <- .C("xDgrsr_arl",as.double(k),
            as.double(g),as.double(zr),as.double(hs),as.double(delta),as.integer(m),
            as.integer(r),as.integer(with0),as.integer(cmode),as.integer(q),
            ans=double(length=1),PACKAGE="spc")$ans
  names(arl) <- "arl"
  return (arl)
}
