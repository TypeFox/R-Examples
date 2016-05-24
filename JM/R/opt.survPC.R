opt.survPC <-
function (thetas) {
    thetas <- relist(thetas, skeleton = list.thetas)
    gammas <- thetas$gammas
    alpha <- thetas$alpha
    Dalpha <- thetas$Dalpha
    xi <- exp(thetas$log.xi)
    eta.tw <- if (!is.null(WW)) as.vector(WW %*% gammas) else 0
    eta.t <- switch(parameterization, 
        "value" = eta.tw + c(WintF.vl %*% alpha) * Y, 
        "slope" = eta.tw + c(WintF.sl %*% Dalpha) * Y.deriv, 
        "both" = eta.tw + c(WintF.vl %*% alpha) * Y + c(WintF.sl %*% Dalpha) * Y.deriv)    
    eta.s <- switch(parameterization, 
        "value" = c(Ws.intF.vl %*% alpha) * Ys, 
        "slope" = c(Ws.intF.sl %*% Dalpha) * Ys.deriv, 
        "both" = c(Ws.intF.vl %*% alpha) * Ys + c(Ws.intF.sl %*% Dalpha) * Ys.deriv)
    log.hazard <- log(xi[ind.D]) + eta.t
    log.survival <- - exp(eta.tw) * rowsum(xi[ind.K] * wkP * exp(eta.s), id.GK, reorder = FALSE)
    dimnames(log.survival) <- NULL
    log.p.tb <- d * log.hazard + log.survival    
    p.bytn <- p.byt * log.p.tb
    -sum(p.bytn %*% wGH, na.rm = TRUE)
}
