##
##  k r i g i n g . R  Kriging Interpolation
##


kriging <- function(u, v, u0, type = c("ordinary", "simple")) {
	stopifnot(is.numeric(u), is.numeric(v), is.numeric(u0))
	if (!is.matrix(u))
		stop("Argument 'u' must be a numeric matrix.")
	n <- nrow(u); m <- ncol(u)
	if (is.vector(v)) {
	    if (length(v) != n)
	    	stop("Length of vector 'v' must be equal to ncol(u).")
	} else if (is.matrix(v)) {
		if (ncol(v) == 1 || nrow(v) == n) {
			v <- c(v)
		} else
		    stop("As matrix 'v' must be a column vector (with ncol(u) elements).")
	} else {
	    stop("Argument 'v' must be a vector or matrix (with ncol(u) elements).")
    }
    if (is.vector(u0)) {
    	if (length(u0) == m) {
    	    u0 <- t(u0)
    	} else {
    	    stop("Length of vector 'u0' must be equal to ncol(u).")
    	}
    } else if (is.matrix(u0)) {
    	if (ncol(u0) != m)
    		stop("Matrix 'u0' must have the same number of colums as 'u'.")
    } else
        stop("Argument 'u0' must be a vector or matrix (with ncol(u) elements).")

    type <- match.arg(type)

    # Define the Variogram
    V  <- distmat(u, u)
    U0 <- distmat(u, u0)

    # Compute kriging formula
    if (type == "simple") {
        w <- v %*% inv(V) %*% U0

    } else if (type == "ordinary") {
    	k <- nrow(u0)
    	C <- matrix(1, n+1, n+1)
    	C[1:n, 1:n] <- V
    	C[n+1, n+1] <- 0

    	D <- matrix(1, n+1, k)
    	D[1:n, 1:k] <- U0

    	v <- c(v, 0)
    	w <- v %*% inv(C) %*% D

    } else  # ntype = 0
        stop("Argument 'type' can only be 'simple' or 'ordinary'.")

    drop(w)     # return as vector
}
