fmodel4pp <-
function(zeta, y, apar, bpar, cpar, dpar, prior = dnorm, ...) {
	if (is.vector(y)) y <- matrix(y, 1, length(y))
	m <- ncol(y)
	prob <- matrix(0, m, 2)
	storage.mode(y) <- "integer"
	storage.mode(prob) <- "double"
  tmp <- .Fortran("fmodel4pp", zeta = as.double(zeta), y = y, 
    m = as.integer(m), s = as.integer(nrow(y)), 
    apar = as.double(apar), bpar = as.double(bpar),
    cpar = as.double(cpar), dpar = as.double(dpar), 
    loglik = as.double(0), prob = prob)
	return(list(post = tmp$loglik + log(prior(zeta, ...)), prob = tmp$prob))
}