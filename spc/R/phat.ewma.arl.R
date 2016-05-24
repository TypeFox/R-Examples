# Computation of EWMA phat ARLs
phat.ewma.arl <- function(lambda, ucl, mu, n, z0, sigma=1, type="known", LSL=-3, USL=3, N=15, qm=25, ntype="coll") {
  if ( lambda <= 0 || lambda > 1 )
    stop("lambda has to be between 0 and 1")
  p.star <- pnorm( LSL ) + pnorm( -USL )
  if ( type == "known" ) {
    if ( ucl <= p.star )
      stop("ucl must be larger than p.star")
  }
  if ( type == "estimated" ) {
    p.star <- 0
    if ( ucl <= 0 )
      stop("ucl must be positive")
  }
  if ( ucl >= 1 )
    stop("ucl must be smaller than 1")
  if ( n < 1 )
    stop("n must be >= 1")
  if ( z0 < p.star | z0 > ucl )
    stop("wrong headstart")
  if ( sigma<1e-12 )
    stop("sigma much too small")
  ctyp <- -1 + pmatch(tolower(type), c("known", "estimated"))
  if ( is.na(ctyp) )
    stop("invalid sigma mode")
  if ( LSL >= USL )
    stop("wrong relationship between lower and upper specification limits (LSL must be smaller than USL)")
  if ( N < 2 )
    stop("N too small")
  if ( qm < 5 )
    stop("qm too small")
  ntyp <- -1 + pmatch(tolower(ntype), c("coll", "markov"))
  if ( is.na(ntyp) )
    stop("wrong label for numerical algorithm")    
  arl <- .C("ewma_phat_arl_coll",
            as.double(lambda), as.double(ucl), as.double(mu), as.double(sigma), as.integer(n),
            as.double(z0), as.integer(ctyp), as.double(LSL), as.double(USL),
            as.integer(N), as.integer(qm), as.integer(ntyp),
            ans=double(length=1), PACKAGE="spc")$ans 
  names(arl) <- "arl"
  arl
}