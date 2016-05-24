##
##  r e p m a t . R
##


repmat <- function(a, n, m = n) {
    if (length(a) == 0) return(c())
    if (!is.numeric(a) && !is.complex(a))
        stop("Argument 'a' must be a numeric or complex.")
    if (is.vector(a))
        a <- matrix(a, nrow = 1, ncol = length(a))
    if (!is.numeric(n) || !is.numeric(m) ||
        length(n) != 1 || length(m) != 1)
        stop("Arguments 'n' and 'm' must be single integers.")
    n <- max(floor(n), 0)
    m <- max(floor(m), 0)
    if (n <= 0 || m <= 0)
        return(matrix(0, nrow = n, ncol = m))

	matrix(1, n, m) %x% a  # Kronecker product
}

Reshape <- function(a, n, m) {
	if (missing(m)) m <- length(a) %/% n
	if (length(a) != n*m)
		stop("Matrix 'a' does not have n*m elements")
	dim(a) <- c(n, m)
	return(a)
}


sortrows <- function(A, k = 1) {  # l <- k+1
    stopifnot(is.numeric(A), is.numeric(k))
    if (!is.matrix(A))
        stop("Argument 'A' must be a numeric matrix.")
    if (length(k) != 1 || floor(k) != ceiling(k))
        stop("Argument 'k' must be a single integer number.")

    n <- nrow(A); m <- ncol(A)
    if (k < 1 || k > m)
        stop("Argument 'k' must satisfy 1 <= k <= ncol(A).")

    o <- order(A[, k])
    A[o, ]
}
