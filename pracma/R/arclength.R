##
##  a r c l e n g t h . R  Arc Length
##


arclength <- function(f, a, b, nmax = 20, tol = 1e-05, ...) {
	stopifnot(is.numeric(a), length(a) == 1,
	          is.numeric(b), length(b) == 1)
 	fun <- match.fun(f)
 	f <- function(x) fun(x, ...)

    fa <- f(a); fb <- f(b)
    m <- length(fa)
    if (length(fa) < 2)
        stop("Argument 'f' must be a parametrized function.")
    if (length(f(c(a, b))) != 2*m)
        stop("Argument 'f' must be a vectorized function.")

    h <- (b - a)
	A <- matrix(0, nmax, nmax)
    A[1, 1] <- sqrt(sum((fb - fa)^2))

    for (i in 1:(nmax-1)) {
    	h <- h/2
    	x <- seq(a, b, by = (b-a)/2^i)
    	y <- c(f(x))
    	X <- matrix(y, ncol = m)
    	dX <- diff(X)

    	A[i+1, 1] <- sum(sqrt(rowSums(dX^2)))
    	for (j in 1:i) {
            A[i+1, j+1] <- (4^j * A[i+1, j] - A[i, j]) / (4^j - 1)
        }
        if (abs(A[i+1, i+1] - A[i, i]) < tol && i > 3) break
    }

    e <- abs(A[i+1, i+1] - A[i, i])
    list(length = A[i+1, i+1], niter = i+1, rel.err = e)
}
