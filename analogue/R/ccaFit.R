rdaFit <- function(X, Y, Z, ...) {
    weight.centre <- function(x, w) {
        w.c <- apply(x, 2, weighted.mean, w = w)
        x <- sweep(x, 2, w.c, "-")
        x <- sweep(x, 1, sqrt(w), "*")
        attr(x, "centre") <- w.c
        x
    }
    ZERO <- 1e-04
    X <- as.matrix(X)
    gran.tot <- sum(X)
    X <- X / gran.tot
    rowsum <- rowSums(X)
    colsum <- colSums(X)
    rc <- outer(rowsum, colsum)
    Xbar <- (X - rc)/sqrt(rc)
    ##tot.chi <- sum(svd(Xbar, nu = 0, nv = 0)$d^2)
    if (!missing(Z) && !is.null(Z)) {
        Z <- as.matrix(Z)
        Z.r <- weight.centre(Z, rowsum)
        Q <- qr(Z.r)
        Z <- qr.fitted(Q, Xbar)
        tmp <- sum(svd(Z, nu = 0, nv = 0)$d^2)
        if (Q$rank) {
            pCCA <- list(rank = Q$rank, tot.chi = tmp, QR = Q,
                         Fit = Z, envcentre = attr(Z.r, "centre"))
            Xbar <- qr.resid(Q, Xbar)
        }
        if (tmp < ZERO)
            pCCA$tot.chi <- 0
    }
}
