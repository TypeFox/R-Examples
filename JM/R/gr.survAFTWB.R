gr.survAFTWB <-
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
    wk.exp.eta.s <- wk * exp(eta.s)
    exp.eta.tw <- exp(eta.tw)
    Vi <- exp.eta.tw * P * rowsum(wk.exp.eta.s, id.GK, reorder = FALSE); dimnames(Vi) <- NULL
    Vii <- d * (sigma.t - 1) / Vi - sigma.t * Vi^(sigma.t - 1)
    scgammas <- - colSums(WW * (d + c((p.byt * Vii * Vi) %*% wGH)), na.rm = TRUE)
    scalpha <- if (parameterization %in% c("value", "both")) {
        rr <- numeric(ncol(WintF.vl))
        for (k in seq_along(rr)) 
            rr[k] <- - sum((p.byt * (d * WintF.vl[, k] * Y + Vii * exp.eta.tw * P * 
                rowsum(wk.exp.eta.s * Ws.intF.vl[, k] * Ys, id.GK, reorder = FALSE))) %*% wGH, na.rm = TRUE)
        rr
    } else NULL
    scalpha.D <- if (parameterization %in% c("slope", "both")) {
        rr <- numeric(ncol(WintF.sl))
        for (k in seq_along(rr)) 
            rr[k] <- - sum((p.byt * (d * WintF.sl[, k] * Y.deriv + Vii * exp.eta.tw * P * 
                rowsum(wk.exp.eta.s * Ws.intF.sl[, k] * Ys.deriv, id.GK, reorder = FALSE))) %*% wGH, na.rm = TRUE)
        rr
    } else NULL
    scsigmat <- if (is.null(scaleWB)) {
         - sigma.t * sum((p.byt * (d / sigma.t + (d - Vi^sigma.t) * log(Vi))) %*% wGH, na.rm = TRUE)
        
    } else NULL
    c(scgammas, scalpha, scalpha.D, scsigmat)
}
