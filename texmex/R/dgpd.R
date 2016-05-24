dgpd <- function(x, sigma, xi, u = 0, log.d=FALSE ) {
  ## record the things below bounds for later
  below.bounds <- x < u
  ## shift and scale
  x <- pmax((x - u) / sigma, 0)

  xix <- .specfun.safe.product(xi, x)
  logrel <- .log1prel(xix) * x

  log.density <- -log(sigma) - log1p(xix) - logrel
  log.density[below.bounds] <- -Inf
  log.density[xix==-1] <- -Inf

  if (!log.d) {
    exp(log.density)
  } else {
    log.density
  }
}


