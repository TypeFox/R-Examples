##
##  l s q l i n . R 
##


lsqlin <- function(A, b, C, d, tol = 1e-13) {
    stopifnot(is.numeric(A), is.numeric(b))
    if (!is.matrix(A))
        stop("Argument 'A' must be a numeric matrix.")

    n <- nrow(A); m <- ncol(A)
    if (is.vector(b)) {
        if (length(b) != n) {
            stop("As vector argument 'b' must have 'nrow(A)' elements.")
        } else {
            b <- as.matrix(b)
            l <- 1
        }
    } else if (is.matrix(b)) {
        if (nrow(b) != n) {
            stop("As Matrix argument 'b' must also have 'nrow(A)' rows.")
        } else {
            l <- ncol(b)
        }
    } else
        stop("Argument 'b' must be a vector or a matrix with n rows.")
    
    if (missing(C) && missing(d)) {
        # x <- pinv(A) %*% b
        # x <- qr.solve(t(A) %*% A, t(A) %*% as.matrix(b))
        x <- qr.solve(A, as.matrix(b))
        return(x)
    } else if ( (missing(C) && !missing(d)) || (!missing(C) && missing(d)))
        stop("Condition 'C * x = d' not fully specified, 'C' or 'd' missing.")

    stopifnot(is.numeric(C), is.numeric(d))
    if (!is.matrix(C) || ncol(C) != m )
        stop("Argument 'C' must be a numeric matrix with 'ncol(C)=ncol(A)'.")

    # xc <- qr.solve(C, d)
    xc <- pinv(C) %*% d
    if ( any(abs(C %*% xc - d) > tol) ) {
        warning("Precondition 'C * x = d' cannot be satisfied (within tolerance 'tol').")
        return(c())
    }

    N <- nullspace(C)                   # (m x k)-matrix, k <= p
    if (is.null(N)) return(c(xc))

    M <- A %*% N                        # (n x k)-matrix
    #xn <- qr.solve(M, b - repmat(A %*% xc, 1, l))
    xn <- pinv(M) %*% (b - repmat(A %*% xc, 1, l))

    x0 <- repmat(xc, 1, l) + N %*% xn
    return(x0)
}
