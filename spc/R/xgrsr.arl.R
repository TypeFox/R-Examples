# Computation of GRSR (Girshick, Rubin, Shiryaev, Roberts) ARLs (mean monitoring)
xgrsr.arl <- function(k, g, mu, zr=0, hs=NULL, sided="one", q=1, MPT=FALSE, r=30) {
  if ( k < 0 )			stop("k has to be non-negative")
  if ( g < 0 )			stop("g has to be positive")
  if ( !is.null(hs) ) {
    if ( hs > g )		stop("wrong headstart")
  } else {
    hs <- 2*g
  }
  q <- round(q)
  if ( q < 1 )			stop("wrong change point position (q)")
  if ( r < 4 )			stop("r is too small")
  ctyp <- pmatch(sided, c("one", "two")) - 1
  if ( is.na(ctyp) )		stop("invalid grsr type")
  arl <- .C("xgrsr_arl",
            as.integer(ctyp), as.double(k), as.double(g),
            as.double(zr), as.double(hs), as.double(mu), as.integer(q), as.integer(r), as.integer(MPT),
            ans=double(length=q), PACKAGE="spc")$ans
  names(arl) <- NULL
  return (arl)
}