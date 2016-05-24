##
##  t r i . R  Triangular matrices
##


tril <- function(M, k = 0) {
    if (k == 0) {
        M[upper.tri(M, diag = FALSE)] <- 0
    } else {
        M[col(M) >= row(M) + k + 1] <- 0
    }
    return(M)
}


triu <- function(M, k = 0) {
    if (k == 0) {
        M[lower.tri(M, diag = FALSE)] <- 0
    } else {
        M[col(M) <= row(M) + k - 1] <- 0
    }
    return(M)
}


##  Format distance matrix
squareform <- function(x) {
    stopifnot(is.numeric(x) || is.complex(x))

    if (is.vector(x)) {
        n <- length(x)
        m <- floor(sqrt(2*n))
        if (m*(m+1) != 2*n)
            stop("Argument 'x' does not correspond to a distance matrix.")
        inds <- c()
        k <- m+1
        for (i in 1:(k-1))
            inds <- c(inds, (1+i+(i-1)*k):(i*k))
        y <- numeric(k*k)
        y[inds] <- x
        y <- matrix(y, k, k) + t(matrix(y, k, k))
    
    } else if (is.matrix(x)) {
        m <- nrow(x); n <- ncol(x)
        if (m != n)
            stop("Argument 'x' must be a vector or a square matrix.")
        if (any(diag(x) != 0))
            stop("Argument 'x' can only have 0s on the diagonal.")
        if (n == 1) return(c())
        inds <- c()
        for (i in 1:(n-1))
            inds <- c(inds, (1+i+(i-1)*n):(i*n))
        y <- x[inds]
    
    } else
        stop("Argument 'x' must be a numeric or complex vector or matrix.")

    return(y)
}
