##
##  p r o c r u s t e s . R  Procrustes Problem
##


procrustes <- function(A, B) {
    stopifnot(is.numeric(A), is.numeric(B))
    if (!is.matrix(A) || !is.matrix(B))
        stop("Arguments 'A' and 'B' must be numeric matrices.")
    if (any(dim(A) != dim(B)))
        stop("Matrices 'A' and 'B' must be of the same size.")

    C <- t(B) %*% A

    Svd <- svd(C)                       # singular value decomposition
    U <- Svd$u
    S <- diag(Svd$d)
    V <- Svd$v

    Q <- U %*% t(V)
    P <- B %*% Q

    R <- A - P;
    r <- sqrt(Trace(t(R) %*% R))        # Frobenius norm: Norm(A - P)
    return(list(P = P, Q = Q, d = r))
}


kabsch <- function(A, B, w = NULL) {
    stopifnot(is.numeric(A), is.numeric(B))
    if ( !is.matrix(A) || !is.matrix(B) )
        stop("Arguments 'A' and 'B' must be numeric matrices.")
    if ( any(dim(A) != dim(B)) )
        stop("Matrices 'A' and 'B' must be of the same size.")

    D <- nrow(A)    # space dimension
    N <- ncol(A)    # number of points
    if (is.null(w)) {
        w <- matrix(1/N, nrow = N, ncol = 1)    # weights as column vector
    } else {
        if (!is.numeric(w) || length(w) != N)
            stop("Argument 'w' must be a (column) vector of length ncol(A).")
    }

    p0 <- A %*% w           # the centroid of A
    q0 <- B %*% w           # the centroid of B
    v1 <- ones(1,N)         # row vector of N ones
    A  <- A - p0 %*% v1     # translating A to center the origin
    B  <- B - q0 %*% v1     # translating B to center the origin

    Pdm <- zeros(D,N)
    for (i in 1:N) Pdm[, i] <- w[i] * A[, i]
    C <- Pdm %*% t(B) 	
    
    Svd <- svd(C)           # singular value decomposition
    V <- Svd$u
    S <- diag(Svd$d)
    W <- Svd$v

    I <- eye(D)
    # more numerically stable than using (det(C) < 0)
    if (det(V %*% t(W)) < 0) I[D, D] <- -1

    U <- W %*% I %*% t(V)
    R <- q0 - U %*% p0

    Diff <- U %*% A - B     # A, B already centered
    lrms <- 0
    for (i in 1:N) lrms <- lrms + w[i] * t(Diff[, i]) %*% Diff[, i]
    lrms <- sqrt(lrms)

    return(list(U = U, R = R, d = lrms))
}
