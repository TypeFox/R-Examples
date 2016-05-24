opt.longPH <-
function (betas) {
    eta.yx <- as.vector(X %*% betas)
    eta.yxT <- as.vector(Xtime %*% betas)
    eta.yxT2 <- as.vector(Xtime2 %*% betas)
    Y <- eta.yxT + Ztime.b
    Y2 <- eta.yxT2 + Ztime2.b
    eta.t <- eta.tw + alpha * Y
    eta.s <- alpha * Y2
    mu.y <- eta.yx + Ztb
    logNorm <- dnorm(y, mu.y, sigma, TRUE)
    log.p.yb <- rowsum(logNorm, id); dimnames(log.p.yb) <- NULL
    log.lambda0T <- log(lambda0[ind.T0])
    log.lambda0T[is.na(log.lambda0T)] <- 0
    log.hazard <- log.lambda0T + eta.t
    S <- matrix(0, n, k)
    S[unq.indT, ] <- rowsum(lambda0[ind.L1] * exp(eta.s), indT, reorder = FALSE)
    log.survival <- - exp(eta.tw) * S
    log.p.tb <- d * log.hazard + log.survival
    p.bytn <- p.byt * (log.p.yb + log.p.tb)
    -sum(p.bytn %*% wGH, na.rm = TRUE)    
}
