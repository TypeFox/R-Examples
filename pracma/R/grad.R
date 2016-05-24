##
##  g r a d . R  Function Gradient
##


grad <- function(f, x0, heps = .Machine$double.eps^(1/3), ...) {
    if (!is.numeric(x0))
        stop("Argument 'x0' must be a numeric value.")

    fun <- match.fun(f)
    f <- function(x) fun(x, ...)

    if (length(f(x0)) != 1)
        stop("Function 'f' must be a univariate function of 2 variables.")
    n <- length(x0)

    hh <- rep(0, n)
    gr <- numeric(n)
    for (i in 1:n) {
        hh[i] <- heps
        gr[i] <- (f(x0 + hh) - f(x0 - hh)) / (2*heps)
        # gr[i] <- (-f(x0+2*hh)+8*f(x0+hh)-8*f(x0-hh)+f(x0-2*hh))/(12*h)
        hh[i] <- 0
    }
    return(gr)
}


# Compute the Jacobian as  J_{ij} = df_i/dx_j  for a vector-valued function
# w/o assuming that the f_i are vectorized.
jacobian <- function(f, x0, heps = .Machine$double.eps^(1/3), ...) {
    if (!is.numeric(x0) || length(x0) == 0)
        stop("Argument 'x' must be a non-empty numeric vector.")

    fun <- match.fun(f)
    f <- function(x) fun(x, ...)

	n <- length(x0)
	m <- length(f(x0))
	jacob <- matrix(NA, m, n)
	hh <- numeric(n)
	for (i in 1:n) {
		hh[i] <- heps
		jacob[, i] <- (f(x0 + hh) - f(x0 - hh)) / (2*heps)
		# jacob[, i] <- (-f(x+2*h)+8*f(x+h)-8*f(x-h)+f(x-2*h))/(12*h.eps)
		hh[i] <- 0
	}
	return(jacob)
}
