`fitPCR` <- function(X, Y, ncomp, n, m) {
    S <- seq_len(ncomp)
    ## model coefficients
    B <- matrix(0, nrow = m, ncol = ncomp)
    Yhat <- matrix(0, nrow = n, ncol = ncomp)
    ## SVD
    SVD <- La.svd(X)
    D <- SVD$d[S]
    TT <- SVD$u[, S, drop = FALSE] %*% diag(D, nrow = ncomp)
    P <- t(SVD$vt[S, , drop = FALSE])
    tQ <- crossprod(TT, Y) / (varExpl <- D^2)
    ## compute coefficients
    for(b in S) {
        bS <- seq_len(b)
        B[, b] <- P[, bS, drop = FALSE] %*% tQ[bS, ]
        Yhat[, b] <- TT[, bS, drop = FALSE] %*% tQ[bS, ]
    }
    list(Yhat = Yhat, B = B, TT = TT, P = P, tQ = tQ, varExpl = varExpl)
}
