##
##  s u b s p a ce . R  Matrix Image


orth <- function(M) {
    if (length(M) == 0)
        return(c())
    if (!is.numeric(M))
        stop("Argument 'M' must be a numeric matrix.")
    if (is.vector(M))
        M <- matrix(c(M), nrow = length(M), ncol = 1)

    svdM <- svd(M)
    U <- svdM$u
    s <- svdM$d
    tol <-  max(dim(M)) * max(s) * .Machine$double.eps

    r <- sum(s > tol)
    U[,1:r, drop = FALSE]
}

subspace <- function(A, B) {
    if (!is.numeric(A) || !is.numeric(B))
        stop("Arguments 'A' and 'B' must be numeric matrices.")
    if (is.vector(A))
        A <- matrix(c(A), nrow = length(A), ncol = 1)
    if (is.vector(B))
        B <- matrix(c(B), nrow = length(B), ncol = 1)
    if (nrow(A) != nrow(B))
        stop("Matrices 'A' and 'B' must have the same number of rows.")

    A <- orth(A)
    B <- orth(B)
    if (ncol(A) < ncol(B)) {
        tmp <- A; A <- B; B <- tmp
    }

    for (k in 1:ncol(A)) {
        B <- B - A[, k] %*% t(A[, k]) %*% B
    }

    asin(min(1, svd(B)$d))
}
