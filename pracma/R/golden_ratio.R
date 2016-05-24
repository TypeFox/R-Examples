##
##  g o l d e n _ r a t i o .R  Golden Ratio Search
##


golden_ratio <- function(f, a, b, ...,
                         maxiter = 100, tol = .Machine$double.eps^0.5)
# Golden Ratio search for a univariate function minimum in a bounded interval
{
    fun <- match.fun(f)
    f <- function(x) fun(x, ...)

	phi <- 1 - (sqrt(5) - 1)/2
	x <- c(a, a + phi*(b-a), b - phi*(b-a), b)
	y2 <- f(x[2])
	y3 <- f(x[3])
	n <- 0
	while (x[3] - x[2] > tol) {
		n <- n + 1
		if (y3 > y2) {
			x <- c(x[1], x[1]+phi*(x[3]-x[1]), x[2], x[3])
			y3 <- y2
			y2 <- f(x[2])
		} else {
			x <- c(x[2], x[3], x[4]-phi*(x[4]-x[2]), x[4])
			y2 <- y3
			y3 <- f(x[3])
		}
		if (n >= maxiter) break
	}
	xm <- (x[2]+x[3])/2
	fxm <- if (abs(f(xm)) <= tol^2) 0.0 else f(xm)
	return(list(xmin=xm, fmin=fxm, iter=n, estim.prec=abs(x[3]-x[2])))
}
