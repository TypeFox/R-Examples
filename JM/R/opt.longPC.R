opt.longPC <-
function (betas) {
    eta.yx <- as.vector(X %*% betas)
    if (parameterization %in% c("value", "both")) {
        Y <- as.vector(Xtime %*% betas) + Ztime.b
        Ys <- as.vector(Xs %*% betas) + Zsb
        eta.t <- eta.tw + c(WintF.vl %*% alpha) * Y
        eta.s <- c(Ws.intF.vl %*% alpha) * Ys
    }
    if (parameterization %in% c("slope", "both")) {
        Y.deriv <- as.vector(Xtime.deriv %*% betas[indFixed]) + Ztime.b.deriv
        Ys.deriv <- as.vector(Xs.deriv %*% betas[indFixed]) + Zsb.deriv
        eta.t <- if (parameterization == "both") 
            eta.t + c(WintF.sl %*% Dalpha) * Y.deriv
        else
            eta.tw + c(WintF.sl %*% Dalpha) * Y.deriv
        eta.s <- if (parameterization == "both")
            eta.s + c(Ws.intF.sl %*% Dalpha) * Ys.deriv 
        else
            c(Ws.intF.sl %*% Dalpha) * Ys.deriv
    }
    mu.y <- eta.yx + Ztb
    logNorm <- dnorm(y, mu.y, sigma, TRUE)
    log.p.yb <- rowsum(logNorm, id)
    log.hazard <- log(xi[ind.D]) + eta.t
    log.survival <- - exp(eta.tw) * rowsum(xi[ind.K] * wkP * exp(eta.s), id.GK, reorder = FALSE)
    dimnames(log.survival) <- NULL
    log.p.tb <- d * log.hazard + log.survival
    p.bytn <- p.byt * (log.p.yb + log.p.tb)
    -sum(p.bytn %*% wGH, na.rm = TRUE)
}
