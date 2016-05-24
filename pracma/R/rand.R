##
##  r a n d . R  Generate Random Matrices
##


rand <- function(n = 1, m = n) {
    stopifnot(is.numeric(n), length(n) <= 2, is.numeric(m))
    if (length(n) == 2)
        return(rand(n[1], n[2]))

    if (length(m) != 1) m <- m[1]
    n <- floor(n)
    m <- floor(m)

    if (n <= 0 || m <= 0) matrix(NA, 0, 0)
    else                  matrix(runif(n*m), nrow=n, ncol=m)
}

randn <- function(n = 1, m = n) {
    stopifnot(is.numeric(n), length(n) <= 2, is.numeric(m))
    if (length(n) == 2)
        return(randn(n[1], n[2]))

    if (length(m) != 1) m <- m[1]
    n <- floor(n)
    m <- floor(m)

    if (n <= 0 || m <= 0) matrix(NA, 0, 0)
    else                  matrix(rnorm(n*m), nrow=n, ncol=m)
}


randi <- function(imax, n = 1, m = n) {
    stopifnot(is.numeric(n), length(n) == 1,
              is.numeric(m), length(m) == 1)
    if (length(imax) == 1) {
        imin <- 1
    } else if (length(imax) == 2) {
        imin <- imax[1]
        imax <- imax[2]
    } else {
        stop("Argument 'imax' must be a scalar or have two elements.")
    }
    if (imin > imax)
        stop("Argument 'imax' must be greater than or equal to 'imin'.")
    n <- floor(n)
    m <- floor(m)

    if (n <= 0 || m <= 0) matrix(NA, 0, 0)
    else matrix(sample(1:imax, n*m, replace=TRUE), nrow=n, ncol=m)
}


rands <- function (n = 1, N = 1, r = 1) 
{
    if (n < 1 || N < 1 || r < 0) return(c())
    X <- randn(n, N+1)
    Y <- sqrt(rowSums(X^2))
    return(r * X/Y)
}


randp <- function(n = 1, r = 1) {
    if (n < 1 || r < 0) return(c())
    x <- rnorm(n); y <- rnorm(n)
    r <- r * sqrt(runif(n)/(x^2 + y^2))
    return(cbind(r*x, r*y))
}


randsample <- function(n, k, w = NULL, replacement = FALSE) {
	stopifnot(is.numeric(n), is.numeric(k))
	if (length(n) == 1) n <- 1:floor(n)
	else                n <- c(n)
	if (k > length(n) && !replacement) {
		warning("k > n or length(n): replacement will be set to TRUE.")
		replacement = TRUE
	}
	if (is.numeric(w)) {
		if (!replacement) replacement = TRUE
		if (length(n) != length(w))
			stop("Weights vector 'w' must have the same length as 'n'.")
	}

	sample(n, k, replace = replacement, prob = w)
}

