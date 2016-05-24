qphat <- function(p, n, mu=0, sigma=1, type="known", LSL=-3, USL=3, nodes=30) {    
  if ( n < 1 )
    stop("n must be >= 1")
  if ( sigma<1e-10 )
    stop("sigma much too small")
  ctyp <- -1 + pmatch(type, c("known", "estimated"))
  if ( is.na(ctyp) )
    stop("invalid sigma mode")
  if ( LSL >= USL )
    stop("wrong relationship between lower and upper specification limits (LSL must be smaller than USL)")
  if ( nodes<2 )
    stop("far too less nodes")

  qf <- rep(NA, length(p))
  for ( i in 1:length(p) ) {
    qf[i] <- NA
    if ( 0<p[i] && p[i]<1 )
      qf[i] <- .C("phat_qf",
                    as.double(p[i]), as.integer(n), as.double(mu), as.double(sigma), as.integer(ctyp),
                    as.double(LSL), as.double(USL), as.integer(nodes),
                    ans=double(length=1), PACKAGE="spc")$ans
  }           
  names(qf) <- "qf"
  qf
}