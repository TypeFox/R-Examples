fmodelrsm <-
function(zeta, y, cpar, dpar, prior = dnorm, ...) {
	if (is.vector(y)) y <- matrix(y, 1, length(y))
  m <- length(dpar)
  r <- length(cpar) + 1
  prob <- matrix(0, m, r)
	storage.mode(y) <- "integer"
  storage.mode(prob) <- "double"
  tmp <- .Fortran("fmodelrsm", zeta = as.double(zeta), y = y,
    m = as.integer(m), r = as.integer(r), s = as.integer(nrow(y)), 
    cpar = as.double(cpar), dpar = as.double(dpar), 
  	loglik = as.double(0), prob = prob)
  return(list(post = tmp$loglik + log(prior(zeta, ...)), prob = tmp$prob))
}
