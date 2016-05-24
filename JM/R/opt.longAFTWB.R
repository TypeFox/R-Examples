opt.longAFTWB <-
function (betas) {
    eta.yx <- as.vector(X %*% betas)
    if (parameterization %in% c("value", "both")) {
        Y <- as.vector(Xtime %*% betas) + Ztime.b
        Ys <- as.vector(Xs %*% betas) + Zsb
        WintF.vl.alph <- c(WintF.vl %*% alpha)
        Ws.intF.vl.alph <- c(Ws.intF.vl %*% alpha)
        eta.t <- eta.tw + WintF.vl.alph * Y
        eta.s <- Ws.intF.vl.alph * Ys
    }
    if (parameterization %in% c("slope", "both")) {
        Y.deriv <- as.vector(Xtime.deriv %*% betas[indFixed]) + Ztime.b.deriv
        Ys.deriv <- as.vector(Xs.deriv %*% betas[indFixed]) + Zsb.deriv
        WintF.sl.alph <- c(WintF.sl %*% Dalpha)
        Ws.intF.sl.alph <- c(Ws.intF.sl %*% Dalpha)
        eta.t <- if (parameterization == "both")
            eta.t + WintF.sl.alph * Y.deriv
        else
            eta.tw + WintF.sl.alph * Y.deriv
        eta.s <- if (parameterization == "both")
            eta.s + Ws.intF.sl.alph * Ys.deriv
        else
            Ws.intF.sl.alph * Ys.deriv
    }
    mu.y <- eta.yx + Ztb
    logNorm <- dnorm(y, mu.y, sigma, TRUE)
    log.p.yb <- rowsum(logNorm, id)
    Vi <- exp(eta.tw) * P * rowsum(wk * exp(eta.s), id.GK, reorder = FALSE); dimnames(Vi) <- NULL
    log.hazard <- log(sigma.t) + (sigma.t - 1) * log(Vi) + eta.t
    log.survival <- - Vi^sigma.t
    log.p.tb <- d * log.hazard + log.survival
    p.bytn <- p.byt * (log.p.yb + log.p.tb)
    -sum(p.bytn %*% wGH, na.rm = TRUE)
}
