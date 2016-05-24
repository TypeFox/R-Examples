LogLik.weibullAFTGH <-
function (thetas) {
    thetas <- relist(thetas, skeleton = list.thetas)
    betas <- thetas$betas
    sigma <- exp(thetas$log.sigma)
    gammas <- thetas$gammas
    alpha <- thetas$alpha
    Dalpha <- thetas$Dalpha
    sigma.t <- if (is.null(scaleWB)) exp(thetas$log.sigma.t) else scaleWB
    D <- thetas$D
    D <- if (diag.D) exp(D) else chol.transf(D)
    eta.yx <- as.vector(X %*% betas)
    eta.tw <- as.vector(WW %*% gammas)
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
    Vi <- exp(eta.tw) * P * rowsum(wk * exp(eta.s), id.GK, reorder = FALSE); dimnames(Vi) <- NULL
    log.hazard <- log(sigma.t) + (sigma.t - 1) * log(Vi) + eta.t
    log.survival <- - Vi^sigma.t
    log.p.tb <- d * log.hazard + log.survival
    log.p.b <- if (control$typeGH == "simple") {
        rep(dmvnorm(b, rep(0, ncz), D, TRUE), each = n)
    } else {
        matrix(dmvnorm(do.call(rbind, lis.b), rep(0, ncz), D, TRUE), n, k, byrow = TRUE)
    }
    p.ytb <- exp(log.p.yb + log.p.tb + log.p.b)
    if (control$typeGH != "simple")
        p.ytb <- p.ytb * VCdets
    dimnames(p.ytb) <- NULL
    p.yt <- c(p.ytb %*% wGH)
    log.p.yt <- log(p.yt)
    - sum(log.p.yt[is.finite(log.p.yt)], na.rm = TRUE)
}
