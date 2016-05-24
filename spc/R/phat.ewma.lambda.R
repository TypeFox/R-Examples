# Computation of EWMA phat lambda minimizing certain out-of-control ARL
phat.ewma.lambda <- function(L0, mu, n, z0, sigma=1, type="known", max_l=1, min_l=.001, LSL=-3, USL=3, qm=25) {
  p.star <- pnorm( LSL/sigma ) + pnorm( -USL/sigma )
  if ( type == "estimated" )
    p.star <- 0
  if ( L0 < 1 )
    stop("L0 is too small")
  if ( n < 1 )
    stop("n must be >= 1")
  if ( z0 < p.star & z0 >= 1 )
    stop("wrong headstart")
  if ( sigma<1e-12 )
    stop("sigma much too small")
  ctyp <- -1 + pmatch(type, c("known", "estimated"))
  if ( is.na(ctyp) )
    stop("invalid sigma mode")
  if ( max_l < min_l | max_l > 1 )
    stop("wrong value for max_l (or min_l)")
  if ( min_l < 1e-4 )
    stop("min_l too small")
  if ( LSL >= USL )
    stop("wrong relationship between lower and upper specification limits (LSL must be smaller than USL)")
  if ( qm < 5 )
    stop("qm too small")
  lambda <- .C("ewma_phat_lambda_coll",
               as.double(L0), as.double(mu), as.double(sigma), as.integer(ctyp),
               as.double(max_l), as.double(min_l), as.integer(n), as.double(z0),
               as.double(LSL), as.double(USL), as.integer(qm),
               ans=double(length=1), PACKAGE="spc")$ans 
  names(lambda) <- "lambda"
  lambda
}