##
##  b r e n t d e k k e r . R  Brent-Dekker Algorithm
##


brentDekker <- function(f, a, b,
                  maxiter = 100, tol = .Machine$double.eps^0.75)
# Brent and Dekker's root finding method,
# based on bisection, secant method and quadratic interpolation
{
    stopifnot(is.numeric(a), is.numeric(b),
              length(a) == 1, length(b) == 1)
    if (!is.function(f) || is.null(f))
        stop("Argument 'f' must be a valid R function.")

	x1 <- a; f1 <- f(x1)
	if (f1 == 0) return(list(root = a, f.root = 0, f.calls = 1, estim.prec = 0))
	x2 <- b; f2 <- f(x2)
	if (f2 == 0) return(list(root = b, f.root = 0, f.calls = 1, estim.prec = 0))
	if (f1*f2 > 0.0)
	    stop("Brent-Dekker: Root is not bracketed in [a, b].")

	x3 <- 0.5*(a+b)
	# Beginning of iterative loop
	niter <- 1
	while (niter <= maxiter) {
		f3 <- f(x3)
		if (abs(f3) < tol) {
		    x0 <- x3
		    break
		}

		# Tighten brackets [a, b] on the root
		if (f1*f3 < 0.0) b <- x3 else a <- x3
		if ( (b-a) < tol*max(abs(b), 1.0) ) {
		    x0 <- 0.5*(a + b)
		    break
	    }

		# Try quadratic interpolation
		denom <- (f2 - f1)*(f3 - f1)*(f2 - f3)
		numer <- x3*(f1 - f2)*(f2 - f3 + f1) + f2*x1*(f2 - f3) + f1*x2*(f3 - f1)
		# if denom==0, push x out of bracket to force bisection
		if (denom == 0) {
			dx <- b - a
		} else {
			dx <- f3*numer/denom
		}

		x <- x3 + dx
		# If interpolation goes out of bracket, use bisection.
		if ((b - x)*(x - a) < 0.0) {
			dx <- 0.5*(b - a)
			x  <- a + dx;
		}

		# Let x3 <-- x & choose new x1, x2 so that x1 < x3 < x2.
		if (x1 < x3) {
			x2 <- x3; f2 <- f3
		} else {
			x1 <- x3; f1 <- f3
		}

		niter <- niter + 1
		if (abs(x - x3) < tol) {
		    x0 <- x
		    break
	    }
		x3 <- x;
	}

    if (niter > maxiter)
        warning("Maximum numer of iterations, 'maxiter', has been reached.")

    prec <- min(abs(x1-x3), abs(x2-x3))
    return(list(root = x0, f.root = f(x0),
                f.calls = niter+2, estim.prec = prec))
}


# alias
brent <- brentDekker

