##
##  c h e b y s h e v . R
##


chebPoly <- function(n, x = NULL) {
    stopifnot(is.numeric(n), length(n) == 1,
              floor(n) == ceiling(n), n >= 0)

	N <- max(2, n+1)
    T <- matrix(0, N, N)

	# Preset degree 0 and 1
	T[1, 1] <- 1
	T[2, 2] <- 1

	# Use recursion formula
	if (n >= 2)
	    for (i in 3:N)
		    T[i, ] <- 2 * c(0, T[i-1, 1:(N-1)]) - T[i-2, ]

    # Reverse with highest coefficient first
    T <- T[ ,ncol(T):1]

    if  (is.null(x)) {
        return(T)
    } else if (is.numeric(x)) {
        return(polyval(T[n+1, ], x))
    } else
        stop("Argument 'x' must be a numeric vector.")
}


chebCoeff <- function(fun, a, b, n) {
    N <- n+1
    # Map interval [a, b] to [-1, 1]
    c1 <- 0.5*(b-a)
    c2 <- 0.5*(b+a)

    # Evaluate function at Chebyshev points (in [a, b])
    k  <- c(0:(N-1))
    y  <- cos(pi*(k+0.5)/N)
    fy <- fun(c1*y + c2)

    # Now compute the Chebyshev coefficients
    c0 <- 2.0 / N
    K  <- matrix(k+0.5, N, 1)
    cheb <- c0 * fy %*% cos(pi/N * K %*% k)

    # Remove too small coefficients and return
    eps <- 2 * .Machine$double.eps    # machine precision
    cheb <- ifelse(abs(cheb) < eps, 0, cheb)
    return(drop(cheb))
}


chebApprox <- function(x, fun, a, b, n) {
	# Compute the Chebyshev polynomials of fun in [a, b]
	cP <- chebPoly(n)
	cC <- chebCoeff(fun, a, b, n)
	p  <- drop(cC %*% cP)
	c0 <- cC[1]

	# Map x into [-1, 1] and evaluate the Chebyshev polynomial
	xx <- (2*x - (b+a))/(b-a)
	yy <- polyval(p, xx) - c0/2
	return(yy)
}
