gr.survPC <-
function (thetas) {
    thetas <- relist(thetas, skeleton = list.thetas)
    gammas <- thetas$gammas
    alpha <- thetas$alpha
    Dalpha <- thetas$Dalpha
    xi <- exp(thetas$log.xi)
    eta.tw <- if (!is.null(WW)) as.vector(WW %*% gammas) else rep(0, n)
    eta.t <- switch(parameterization, 
        "value" = eta.tw + c(WintF.vl %*% alpha) * Y, 
        "slope" = eta.tw + c(WintF.sl %*% Dalpha) * Y.deriv, 
        "both" = eta.tw + c(WintF.vl %*% alpha) * Y + c(WintF.sl %*% Dalpha) * Y.deriv)    
    exp.eta.s <- exp(switch(parameterization, 
        "value" = c(Ws.intF.vl %*% alpha) * Ys, 
        "slope" = c(Ws.intF.sl %*% Dalpha) * Ys.deriv, 
        "both" = c(Ws.intF.vl %*% alpha) * Ys + c(Ws.intF.sl %*% Dalpha) * Ys.deriv))
    exp.eta.tw <- exp(eta.tw)
    Int <-  wkP * exp.eta.s
    Int2 <- xi[ind.K] * Int
    scgammas <- if (!is.null(WW)) {
        - colSums(WW * (d - c((p.byt * (exp.eta.tw * 
            rowsum(Int2, id.GK, reorder = FALSE))) %*% wGH)), na.rm = TRUE)
    } else NULL
    scalpha <- if (parameterization %in% c("value", "both")) {
        rr <- numeric(ncol(WintF.vl))
        for (k in seq_along(rr)) 
            rr[k] <- - sum((p.byt * (d * WintF.vl[, k] * Y - exp.eta.tw * 
                rowsum(Int2 * Ws.intF.vl[, k] * Ys, id.GK, reorder = FALSE))) %*% wGH, na.rm = TRUE)
        rr
    } else NULL
    scalpha.D <- if (parameterization %in% c("slope", "both")) {
        rr <- numeric(ncol(WintF.sl))
        for (k in seq_along(rr)) 
            rr[k] <- - sum((p.byt * (d * WintF.sl[, k] * Y.deriv - exp.eta.tw * 
                rowsum(Int2 * Ws.intF.sl[, k] * Ys.deriv, id.GK, reorder = FALSE))) %*% wGH, na.rm = TRUE)
        rr
    } else NULL
    scxi <- numeric(Q)
    for (i in 1:Q) {
        i1 <- ind.D == i
        i2 <- ind.K == i
        i3 <- ind.D >= i
        ki <- c((p.byt[i3, ] * (exp.eta.tw[i3] * rowsum(Int[i2, ], id.GK[i2], reorder = FALSE))) %*% wGH)
        kk <- numeric(n); kk[i3] <- ki
        scxi[i] <- - xi[i] * sum((d * i1)/xi[i] - kk)
    }
    c(scgammas, scalpha, scalpha.D, scxi)
}
