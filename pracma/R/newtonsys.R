##
##  n e w t o n s y s . R
##


# Finds a zero of a nonlinear system by Newton's method
newtonsys <- function(Ffun, x0, Jfun = NULL, ...,
	                  maxiter = 100, tol = .Machine$double.eps^(1/2)) {
	.vnorm <- function(x) { sqrt(sum(x^2)) }
	if (is.null(Jfun)) Jfun <- function(x, ...) jacobian(Ffun, x, ...)

	niter <- 0; err <- tol + 1
	x <- x0
	while (err >= tol && niter < maxiter) {
		niter <- niter + 1
		F <- Ffun(x, ...)
		J <- Jfun(x, ...)
		delta <- -1 * solve(J, F)
		x <- x + delta
		err <- .vnorm(delta)
	}
	F <- .vnorm(Ffun(x, ...))
	if (niter > maxiter && err > tol) {
		cat("Fails to converge within maximum number of iterations.\n",
		    "The iterate returned has relative residual ", F, "\n", sep="")
	}
	return(list(zero=x, fnorm=F, niter=niter))
}
