fmodelnrm <-
function(zeta, y, apar, bpar, prior = dnorm, ...) {
	if (is.vector(y)) y <- matrix(y, 1, length(y))
	m <- ncol(y)
  r <- ncol(apar)
  prob <- matrix(0, m, r)
	storage.mode(y) <- "integer"
  storage.mode(apar) <- "double"
  storage.mode(bpar) <- "double"
  storage.mode(prob) <- "double"
  tmp <- .Fortran("fmodelnrm", zeta = as.double(zeta), y = y,
    m = as.integer(m), r = as.integer(r), s = as.integer(nrow(y)),
  	apar = apar, bpar = bpar, 
    loglik = as.double(0), prob = prob)
  return(list(post = tmp$loglik + log(prior(zeta, ...)), prob = tmp$prob))
}