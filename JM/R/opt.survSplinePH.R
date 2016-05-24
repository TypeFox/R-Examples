opt.survSplinePH <-
function (thetas) {
    thetas <- relist(thetas, skeleton = list.thetas)
    gammas <- thetas$gammas
    alpha <- thetas$alpha
    Dalpha <- thetas$Dalpha
    gammas.bs <- thetas$gammas.bs
    eta.tw1 <- if (!is.null(W1)) as.vector(W1 %*% gammas) else rep(0, n)
    eta.tw2 <- as.vector(W2 %*% gammas.bs)
    eta.t <- switch(parameterization, 
        "value" = eta.tw2 + eta.tw1 + c(WintF.vl %*% alpha) * Y, 
        "slope" = eta.tw2 + eta.tw1 + c(WintF.sl %*% Dalpha) * Y.deriv, 
        "both" = eta.tw2 + eta.tw1 + c(WintF.vl %*% alpha) * Y + 
            c(WintF.sl %*% Dalpha) * Y.deriv)    
    eta.s <- switch(parameterization, 
        "value" = c(Ws.intF.vl %*% alpha) * Ys,
        "slope" = c(Ws.intF.sl %*% Dalpha) * Ys.deriv, 
        "both" = c(Ws.intF.vl %*% alpha) * Ys + 
            c(Ws.intF.sl %*% Dalpha) * Ys.deriv)
    eta.ws <- as.vector(W2s %*% gammas.bs)
    log.hazard <- eta.t
    log.survival <- - exp(eta.tw1) * P * rowsum(wk * exp(eta.ws + eta.s), 
        id.GK, reorder = FALSE)
    dimnames(log.survival) <- NULL
    log.p.tb <- rowsum(d * log.hazard + log.survival, idT, reorder = FALSE)
    p.bytn <- p.byt * log.p.tb
    -sum(p.bytn %*% wGH, na.rm = TRUE)
}
