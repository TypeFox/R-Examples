opt.survWB <-
function (thetas) {
    thetas <- relist(thetas, skeleton = list.thetas)
    gammas <- thetas$gammas
    alpha <- thetas$alpha
    Dalpha <- thetas$Dalpha
    sigma.t <- if (is.null(scaleWB)) exp(thetas$log.sigma.t) else scaleWB
    eta.tw <- as.vector(WW %*% gammas) 
    eta.t <- switch(parameterization, 
        "value" = eta.tw + c(WintF.vl %*% alpha) * Y, 
        "slope" = eta.tw + c(WintF.sl %*% Dalpha) * Y.deriv,
        "both" = eta.tw + c(WintF.vl %*% alpha) * Y + c(WintF.sl %*% Dalpha) * Y.deriv)    
    eta.s <- switch(parameterization, 
        "value" = c(Ws.intF.vl %*% alpha) * Ys,
        "slope" = c(Ws.intF.sl %*% Dalpha) * Ys.deriv, 
        "both" = c(Ws.intF.vl %*% alpha) * Ys + c(Ws.intF.sl %*% Dalpha) * Ys.deriv)
    log.hazard <- log(sigma.t) + (sigma.t - 1) * logT + eta.t
    log.survival <- - exp(eta.tw) * P * rowsum(wk * exp(log(sigma.t) + 
        (sigma.t - 1) * log.st + eta.s), id.GK, reorder = FALSE)
    dimnames(log.survival) <- NULL
    log.p.tb <- d * log.hazard + log.survival    
    p.bytn <- p.byt * log.p.tb
    -sum(p.bytn %*% wGH, na.rm = TRUE)
}
