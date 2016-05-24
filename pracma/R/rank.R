##
##  r a n k . R  Matrix Rank
##


Rank <- function(M) {
    if (length(M) == 0)
        return(0)
    if (!is.numeric(M))
        stop("Argument 'M' must be a numeric matrix.")
    if (is.vector(M))
        M <- matrix(c(M), nrow = length(M), ncol = 1)

    # The MASS way
    r1 <- qr(M)$rank

    # The Matlab way
    sigma <- svd(M)$d
    tol <- max(dim(M)) * max(sigma) * .Machine$double.eps
    r2 <- sum(sigma > tol)

    if (r1 != r2)
        warning("Rank calculation may be problematic.")
    return(r2)
}


nullspace <- function(M) {
    if (!is.numeric(M))
        stop("Argument 'M' must be a numeric matrix.")
    if (is.vector(M))
        M <- matrix(c(M), nrow = length(M), ncol = 1)

    qrM <- qr(t(M))
    rnk <- qrM$rank
    if (rnk == ncol(M)) return(NULL)

    inds <- if (rnk == 0) 1:ncol(M) else -(1:rnk)
    qrQ <- qr.Q(qrM, complete = TRUE)[, inds, drop = FALSE]

    if (length(qrQ) == 0) return(NULL)
    else                  return(qrQ)
}


null <- nullspace
