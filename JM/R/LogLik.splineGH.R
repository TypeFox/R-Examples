LogLik.splineGH <-
function (thetas) {
    thetas <- relist(thetas, skeleton = list.thetas)
    betas <- thetas$betas
    sigma <- exp(thetas$log.sigma)
    gammas <- thetas$gammas
    gammas.bs <- thetas$gammas.bs
    alpha <- thetas$alpha
    Dalpha <- thetas$Dalpha
    D <- thetas$D
    D <- if (diag.D) exp(D) else chol.transf(D)
    eta.yx <- as.vector(X %*% betas)
    eta.tw1 <- if (!is.null(W1)) as.vector(W1 %*% gammas) else rep(0, n)
    eta.tw2 <- as.vector(W2 %*% gammas.bs)
    if (parameterization %in% c("value", "both")) {
        Y <- as.vector(Xtime %*% betas) + Ztime.b
        Ys <- as.vector(Xs %*% betas) + Zsb
        eta.t <- eta.tw2 + eta.tw1 + c(WintF.vl %*% alpha) * Y
        eta.s <- c(Ws.intF.vl %*% alpha) * Ys
    }
    if (parameterization %in% c("slope", "both")) {
        Y.deriv <- as.vector(Xtime.deriv %*% betas[indFixed]) + Ztime.b.deriv
        Ys.deriv <- as.vector(Xs.deriv %*% betas[indFixed]) + Zsb.deriv
        eta.t <- if (parameterization == "both")
            eta.t + c(WintF.sl %*% Dalpha) * Y.deriv
        else
            eta.tw2 + eta.tw1 + c(WintF.sl %*% Dalpha) * Y.deriv
        eta.s <- if (parameterization == "both")
            eta.s + c(Ws.intF.sl %*% Dalpha) * Ys.deriv
        else
            c(Ws.intF.sl %*% Dalpha) * Ys.deriv
    }
    eta.ws <- as.vector(W2s %*% gammas.bs)
    mu.y <- eta.yx + Ztb
    logNorm <- dnorm(y, mu.y, sigma, TRUE)
    log.p.yb <- rowsum(logNorm, id)    
    log.hazard <- eta.t
    log.survival <- - exp(eta.tw1) * P * rowsum(wk * exp(eta.ws + eta.s), 
        id.GK, reorder = FALSE)
    dimnames(log.survival) <- NULL
    log.p.tb <- rowsum(d * log.hazard + log.survival, idT, reorder = FALSE)
    log.p.b <- if (control$typeGH == "simple") {
        rep(dmvnorm(b, rep(0, ncz), D, TRUE), each = n)
    } else {
        matrix(dmvnorm(do.call(rbind, lis.b), rep(0, ncz), D, TRUE), 
            n, k, byrow = TRUE)
    }
    p.ytb <- exp(log.p.yb + log.p.tb + log.p.b)
    if (control$typeGH != "simple")
        p.ytb <- p.ytb * VCdets
    dimnames(p.ytb) <- NULL
    p.yt <- c(p.ytb %*% wGH)
    log.p.yt <- log(p.yt)
    - sum(log.p.yt[is.finite(log.p.yt)], na.rm = TRUE)
}
