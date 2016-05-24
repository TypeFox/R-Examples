##
##  f i b s e a r c h . R  Fibonacci Search
##


fibsearch <- function(f, a, b, ...,
					  endp = FALSE, tol = .Machine$double.eps^(1/2))
# Fibonacci search for a univariate function minimum in a bounded interval
{
	if (a >= b) stop("Left endpoint a must be smaller than b.")
	tol <- max(tol, .Machine$double.eps)

	# Compute Fibonacci numbers [F0,] F1, F2, ..., Fm
	F <- c(1, 2); n <- 2
	while (F[n] <= 2*(b-a)/tol) { F[n+1] <- F[n] + F[n-1]; n <- n + 1 }

	# Initialize values (k == 0)
	x1 <- a; x2 <- b
	xa <- a + (b-a) * F[n-2]/F[n]; fxa <- f(xa, ...)
	xb <- a + (b-a) * F[n-1]/F[n]; fxb <- f(xb, ...)

	# Compute iteration
	k <- 1
	while (k <= n-3 && xa < xb && (x2 - x1) >= tol) {
		if (fxa > fxb) {
			x1 <- xa; xa <- xb
			xb <- x1 + (x2-x1) * F[n-k-1]/F[n-k]
			fxa <- fxb; fxb <- f(xb, ...)
		} else {
			x2 <- xb; xb <- xa
			xa <- x1 + (x2-x1) * F[n-k-2]/F[n-k]
			fxb <- fxa; fxa <- f(xa, ...)
		}
		k <- k + 1
	}

	# Finally use the mean and consider endpoints
	xmin <- (xa+xb)/2; fmin <- f(xmin, ...)
	if (endp) {
		fa <- f(a, ...); fb <- f(b, ...)
		if (xmin-a < tol && fa < fmin) {
			xmin <- a; fmin <- fa
		} else {
			if (b - xmin < tol && fb < fmin) {
				xmin <- b; fmin <- fb
			}
		}
	}
	estim.prec <- max(xmin-x1, x2-xmin)
	return(list(xmin=xmin, fmin=fmin, niter=k, estim.prec=estim.prec))
}
