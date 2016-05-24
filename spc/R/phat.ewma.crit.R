# Computation of EWMA phat upper control limits
phat.ewma.crit <- function(lambda, L0, mu, n, z0, sigma=1, type="known", LSL=-3, USL=3, N=15, qm=25) {
  if ( lambda <= 0 || lambda > 1 )
    stop("lambda has to be between 0 and 1")
  p.star <- pnorm( LSL/sigma ) + pnorm( -USL/sigma )
  if ( type == "estimated" )
    p.star <- 0
  if ( L0 < 1 )
    stop("L0 is too small")
  if ( n < 1 )
    stop("n must be >= 1")
  if ( z0 < p.star & z0 >= 1 )
    stop("wrong headstart")    
  if ( sigma<1e-10 )
    stop("sigma much too small")
  ctyp <- -1 + pmatch(type, c("known", "estimated"))
  if ( is.na(ctyp) )
    stop("invalid sigma mode")
  if ( LSL >= USL )
    stop("wrong relationship between lower and upper specification limits (LSL must be smaller than USL)")
  if ( N < 3 )
    stop("N too small")
  if ( qm < 5 )
    stop("qm too small")
  ucl <- .C("ewma_phat_crit_coll",
            as.double(lambda), as.double(L0), as.double(mu), as.double(sigma), as.integer(n),
            as.double(z0), as.integer(ctyp), as.double(LSL), as.double(USL), as.integer(N), as.integer(qm),
            ans=double(length=1), PACKAGE="spc")$ans 
  names(ucl) <- "ucl"
  ucl
}