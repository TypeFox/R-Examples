dgev <- function(x, mu, sigma, xi, log.d=FALSE){
  ## shift and scale
  x <- (x - mu) / sigma

  xix <- .specfun.safe.product(xi, x)
  logrel <- .log1prel(xix) * x

  log.density <- -log(sigma) - log1p(xix) - logrel - exp(-logrel)
  ## make exp(Inf) > Inf
  log.density[logrel==(-Inf)] <- -Inf

  if (!log.d) {
    exp(log.density)
  } else {
    log.density
  }
}
