gr.survSplinePH <-
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
    exp.eta.tw.P <- exp(eta.tw1) * P
    Int <- wk * exp(eta.ws + eta.s)
    scgammas1 <- if (!is.null(W1)) {
        ki <- exp.eta.tw.P * rowsum(Int, id.GK, reorder = FALSE)
        scg1 <- numeric(ncol(W1))
        for (jj in seq_along(scg1)) {
            tt <- rowsum(W1[, jj] * ki, idT, reorder = FALSE)
            scg1[jj] <- sum(c((p.byt * tt) %*% wGH), na.rm = TRUE)
        }
        - colSums(W1 * d, na.rm = TRUE) + scg1
    } else 
        NULL
    scgammas2 <- numeric(nk)
    for (i in 1:nk) {
        kk <- exp.eta.tw.P * rowsum(Int * W2s[, i], id.GK, reorder = FALSE)
        kk <- rowsum(kk, idT, reorder = FALSE)
        scgammas2[i] <- - sum(W2[, i] * d) + sum(c((p.byt * kk) %*% wGH))
    }
    scalpha <- if (parameterization %in% c("value", "both")) {
        rr <- numeric(ncol(WintF.vl))
        for (k in seq_along(rr)) {
            rrr <- exp.eta.tw.P * rowsum(Int * Ws.intF.vl[, k] * Ys, 
                    id.GK, reorder = FALSE)
            rrr <- rowsum(rrr, idT, reorder = FALSE)            
            rr[k] <- - sum((p.byt * (rowsum(d * WintF.vl[, k] * Y, idT, 
                    reorder = FALSE) - rrr)) %*% wGH, na.rm = TRUE)
        }
        rr
    } else NULL
    scalpha.D <- if (parameterization %in% c("slope", "both")) {
        rr <- numeric(ncol(WintF.sl))
        for (k in seq_along(rr)) {
            rrr <- exp.eta.tw.P * rowsum(Int * Ws.intF.sl[, k] * 
                Ys.deriv, id.GK, reorder = FALSE)
            rrr <- rowsum(rrr, idT, reorder = FALSE)            
            rr[k] <- - sum((p.byt * (rowsum(d * WintF.sl[, k] * Y.deriv, 
                idT, reorder = FALSE) - rrr)) %*% wGH, na.rm = TRUE)
        }
        rr
    } else NULL
    c(scgammas1, scalpha, scalpha.D, scgammas2)
}
