##
##  q r . R  QR Factorization
##


# Modified Gram-Schmidt process
gramSchmidt <- function(A, tol = .Machine$double.eps^0.5) {
    stopifnot(is.numeric(A), is.matrix(A))
    m <- nrow(A); n <- ncol(A)
    if (m < n)
        stop("No. of rows of 'A' must be greater or equal no. of colums.")

    Q <- matrix(0, m, n)
    R <- matrix(0, n, n)
    for (k in 1:n) {
        Q[, k] <- A[, k]
        if (k > 1) {
            for (i in 1:(k-1)) {
                R[i, k] <- t(Q[, i]) %*% Q[, k]
                Q[ , k] <- Q[, k] - R[i, k] * Q[, i]
            }
        }
        R[k, k] <- Norm(Q[, k])
        if (abs(R[k, k]) <= tol)
            stop("Matrix 'A' does not have full rank.")
        Q[, k] <- Q[, k] / R[k, k]
    }
    return(list(Q = Q, R = R))
}

qrSolve <- function(A, b) {
    stopifnot(is.numeric(A), is.matrix(A), is.numeric(b))
    m <- nrow(A); n <- ncol(A)
    b <- c(b)
    if (m < n || length(b) != m)
        stop("Matrix 'A' and vektor 'b' have non-fitting dimensions.")

    gs <- householder(A)
    Q <- gs$Q; R <- gs$R

    b <- t(Q[, 1:n]) %*% b
    x <- numeric(n)
    x[n] <- b[n] / R[n, n]
    for (k in (n-1):1)
        x[k] <- (b[k] - R[k, (k+1):n] %*% x[(k+1):n]) / R[k, k]
    return(x)
}


# Givens transformation
.givens <- function(xk, xl) {
    if (xl != 0) {
        r <- Norm(c(xk, xl))
        G <- matrix(c(xk, -xl, xl, xk), 2, 2) / r
        x <- as.matrix(c(r, 0))
    } else {
        G <- eye(2)
        x <- as.matrix(c(xk, 0))
    }
    return(list(G = G, x = x))
}

# Givens QR decomposition
givens <- function(A) {  # n >= m
    stopifnot(is.numeric(A), is.matrix(A))
    n <- nrow(A); m <- ncol(A)
    if (n != m)
        stop("Matrix 'A' must be a square matrix.")

    Q <- eye(n)
    for (k in 1:(n-1)) {
        l <- which.max(abs(A[(k+1):n, k])) + k
        if (A[k, k] == 0 && A[l, k] == 0)
            stop("Matrix 'A' does not have full rank.")
        j <- which(A[(k+1):n, k] != 0) + k
        j <- unique(c(l, j[j != 1]))
        for (h in j) {
            gv <- .givens(A[k, k], A[h, k])
            G <- gv$G; x <- gv$x
            Q[c(k, h), ] <- G %*% Q[c(k, h), ]
            A[k, k] <- x[1]
            A[h, k] <- 0
            A[c(k, h), (k+1):m] <- G %*% A[c(k, h), (k+1):m]
        }
    }
    return(list(Q = t(Q), R = triu(A)))
}


# Householder transformation
householder <- function(A) {
    m <- nrow(A); n <- ncol(A)
    Q <- eye(m)
    for (k in 1:min(m-1, n)) {
        ak <- A[k:m, k, drop = FALSE]
        s  <- if (ak[1] >= 0) 1 else -1
        vk <- ak + s * Norm(ak) * c(1, rep(0, m-k))
        vk2 <- c(t(vk) %*% vk)
        Hk <- eye(m-k+1) - 2/vk2 * (vk %*% t(vk))
        if (k == 1) Qk <- Hk
        else        Qk <- blkdiag(eye(k-1), Hk)
        A <- Qk %*% A
        Q <- Q %*% Qk
    }
    return(list(Q = Q, R = A))
}
