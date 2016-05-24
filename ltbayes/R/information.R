information <- function(fmodel, y, zeta, observed = FALSE, ...) {
	if (!is.vector(y) & !is.matrix(y)) stop("y must be a vector or matrix of item responses.")
	if (is.vector(y)) y <- matrix(y, 1, length(y))
	if (observed) {
		if (!missing(zeta)) warning("Observed information computed at MLE only; zeta argument ignored.")
		zeta.mode <- postmode(fmodel, y, prior = function(z) return(1), ...)$zeta
		inf <- -numDeriv::hessian(function(z, ...) fmodel(z, ...)$post, zeta.mode, y = y,
			prior = function(z) return(1), ...)
		return(list(test = inf))
	}
	else {
		if (missing(zeta)) zeta <- postmode(fmodel, y, prior = function(z) return(1), ...)$zeta
		if (nrow(y) > 1) stop("Fisher information only available for single response patterns.")
		p <- fmodel(zeta, y, ...)$prob
		m <- nrow(p)
		r <- ncol(p)
		g <- matrix(NA, m, r)
		h <- matrix(NA, m, r)
		for (j in 1:m) {
			for (k in 1:r) {
				f <- function(z, ...) fmodel(z, ...)$prob[j,k]
				g[j,k] <- numDeriv::grad(f, zeta, y = y, ...)
				h[j,k] <- numDeriv::hessian(f, zeta, y = y, ...)
			}
		}
		inf <- g^2/p - h
		return(list(test = sum(inf), item = apply(inf, 1, sum), category = inf))
	}
}