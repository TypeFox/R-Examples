pphat <- function(q, n, mu=0, sigma=1, type="known", LSL=-3, USL=3, nodes=30) {    
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
    
  p.star <- pnorm( LSL/sigma ) + pnorm( -USL/sigma )
  if ( type == "estimated" ) p.star <- 0

  cdf <- rep(NA, length(q))
  for ( i in 1:length(q) ) {
    cdf[i] <- 0
    if ( q[i] >= 1 ) cdf[i] <- 1
    if ( p.star<q[i] && q[i]<1 )
      cdf[i] <- .C("phat_cdf",
                    as.double(q[i]), as.integer(n), as.double(mu), as.double(sigma), as.integer(ctyp),
                    as.double(LSL), as.double(USL), as.integer(nodes),
                    ans=double(length=1), PACKAGE="spc")$ans
  }           
  names(cdf) <- NULL
  cdf
}