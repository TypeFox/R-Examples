gr.longPH <-
function (betas) {
    eta.yx <- as.vector(X %*% betas)
    eta.yxT <- as.vector(Xtime %*% betas)
    eta.yxT2 <- as.vector(Xtime2 %*% betas)
    Y <- eta.yxT + Ztime.b
    Y2 <- eta.yxT2 + Ztime2.b
    eta.t <- eta.tw + alpha * Y
    eta.s <- alpha * Y2
    exp.eta.s <- exp(eta.s)
    exp.eta.tw <- exp(eta.tw)
    sc1 <- - crossprod(X, y - eta.yx - Zb) / sigma^2
    Int <- lambda0[ind.L1] * exp(eta.s) * alpha
    sc2 <- numeric(ncx)
    for (i in 1:ncx) {
        S <- matrix(0, n, k)
        S[unq.indT, ] <- rowsum(Int * Xtime2[, i], indT, reorder = FALSE)
        ki <- exp.eta.tw * S
        kii <- c((p.byt * ki) %*% wGH)
        sc2[i] <- - sum(d * alpha * Xtime[, i] - kii, na.rm = TRUE)
    }
    c(sc1 + sc2)
}
