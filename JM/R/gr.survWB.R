gr.survWB <-
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
    exp.eta.tw.P <- exp(eta.tw) * P
    Int <- wk * exp(log(sigma.t) + (sigma.t - 1) * log.st + eta.s)
    ki <- exp.eta.tw.P * rowsum(Int, id.GK, reorder = FALSE)
    kii <- c((p.byt * ki) %*% wGH)
    scgammas <- - colSums(WW * (d - kii), na.rm = TRUE)
    scalpha <- if (parameterization %in% c("value", "both")) {
        rr <- numeric(ncol(WintF.vl))
        for (k in seq_along(rr)) 
            rr[k] <- - sum((p.byt * (d * WintF.vl[, k] * Y - exp.eta.tw.P * 
                rowsum(Int * Ws.intF.vl[, k] * Ys, id.GK, reorder = FALSE))) %*% wGH, na.rm = TRUE)
        rr
    } else NULL
    scalpha.D <- if (parameterization %in% c("slope", "both")) {
        rr <- numeric(ncol(WintF.sl))
        for (k in seq_along(rr)) 
            rr[k] <- - sum((p.byt * (d * WintF.sl[, k] * Y.deriv - exp.eta.tw.P * 
                rowsum(Int * Ws.intF.sl[, k] * Ys.deriv, id.GK, reorder = FALSE))) %*% wGH, na.rm = TRUE)
        rr
    } else NULL
    scsigmat <- if (is.null(scaleWB)) {
        Int2 <- st^(sigma.t - 1) * (1 + sigma.t * log.st) * exp(eta.s)
         - sigma.t * sum((p.byt * (d * (1/sigma.t + logT) - exp.eta.tw.P * 
            rowsum(wk * Int2, id.GK, reorder = FALSE))) %*% wGH, na.rm = TRUE)
    } else NULL
    c(scgammas, scalpha, scalpha.D, scsigmat)
}
