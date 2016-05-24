`fitCoinertia` <- function(X, Dp, Y, Dq, Dn, n.axes) {
    ax.names <- paste("COCA", 1:n.axes, sep = " ")
    Dp05 <- diag(sqrt(Dp))
    Dq05 <- diag(sqrt(Dq))
    A <- Dp05 %*% t(X) %*% diag(Dn) %*% Y %*% Dq05
    svdA <- La.svd(A)
    U <- diag(1 / sqrt(Dp)) %*% svdA$u
    V <- diag(1 / sqrt(Dq)) %*% t(svdA$vt)
    Ksi <- X %*% diag(Dp) %*% U
    Psi <- Y %*% diag(Dq) %*% V
    L <- diag(svdA$d)
    L <- L * L
    seqA <- seq_len(n.axes)
    U1 <- U[, seqA, drop = FALSE]
    U2 <- V[, seqA, drop = FALSE]
    colnames(U1) <- colnames(U2) <- ax.names
    rownames(U1) <- colnames(X)
    rownames(U2) <- colnames(Y)
    X1 <- Ksi[, seqA, drop = FALSE]
    X2 <- Psi[, seqA, drop = FALSE]
    colnames(X1) <- colnames(X2) <- ax.names
    lambda <- diag(L[seqA, seqA, drop = FALSE])
    names(lambda) <- ax.names
    retval <- list(scores = list(species = list(Y = U1, X= U2),
                   site = list(Y = X1, X = X2)),
                   lambda = lambda, n.axes = n.axes, call = match.call())
    class(retval) <- "fitCoinertia"
    retval
}

