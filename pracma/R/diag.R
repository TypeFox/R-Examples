##
##  d i a g . R  Matrix diagonal
##


Diag <- function(x, k=0) {
    if (!is.numeric(x) && !is.complex(x))
        stop("Argument 'x' must be a real or complex vector or matrix.")
    if (!is.numeric(k) || k != round(k))
        stop("Argument 'k' must be an integer.")

    # if (length(x) == 1) return(x)
    if (is.matrix(x)) {
        n <- nrow(x); m <- ncol(x)
        if (k >= m || -k >= n) {
            y <- matrix(0, nrow=0, ncol=0)
        } else {
            y <- x[col(x) == row(x) + k]
        }
    } else {
        if (is.vector(x)) {
            n <- length(x)
            m <- n + abs(k)
            y <- matrix(0, nrow=m, ncol=m)
            y[col(y) == row(y) + k] <- x
        } else {
            stop("Argument 'x' must be a real or complex vector or matrix.")
        }
    }
    return(y)
}
