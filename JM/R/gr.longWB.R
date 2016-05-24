gr.longWB <-
function (betas) {
    eta.yx <- as.vector(X %*% betas)
    if (parameterization %in% c("value", "both")) {
        Ys <- as.vector(Xs %*% betas) + Zsb
        WintF.vl.alph <- c(WintF.vl %*% alpha)
        Ws.intF.vl.alph <- c(Ws.intF.vl %*% alpha)
        eta.s <- Ws.intF.vl.alph * Ys
    }
    if (parameterization %in% c("slope", "both")) {
        Ys.deriv <- as.vector(Xs.deriv %*% betas[indFixed]) + Zsb.deriv
        WintF.sl.alph <- c(WintF.sl %*% Dalpha)
        Ws.intF.sl.alph <- c(Ws.intF.sl %*% Dalpha)
        eta.s <- if (parameterization == "both")
            eta.s + Ws.intF.sl.alph * Ys.deriv
        else
            Ws.intF.sl.alph * Ys.deriv
    }
    exp.eta.tw.P <- exp(eta.tw) * P
    sc1 <- - crossprod(X, y - eta.yx - Zb) / sigma^2
    Int <- wk * exp(log(sigma.t) + (sigma.t - 1) * log.st + eta.s)
    sc2 <- numeric(ncx)
    for (i in 1:ncx) {
        ki <- exp.eta.tw.P * switch(parameterization,
            "value" = rowsum(Int * Ws.intF.vl.alph * Xs[, i], id.GK, reorder = FALSE),
            "slope" = {ii <- match(i, indFixed); 
                if (is.na(ii)) 0 else rowsum(Int * Ws.intF.sl.alph * Xs.deriv[, ii], id.GK, reorder = FALSE)},
            "both" = {ii <- match(i, indFixed); 
                rowsum(Int * (Ws.intF.vl.alph * Xs[, i] + 
                    Ws.intF.sl.alph * if (is.na(ii)) 0 else Xs.deriv[, ii]), id.GK, reorder = FALSE)} 
        )
        kii <- c((p.byt * ki) %*% wGH)
        sc2[i] <- switch(parameterization,
            "value" = - sum(d * WintF.vl.alph * Xtime[, i] - kii, na.rm = TRUE),
            "slope" = {ii <- match(i, indFixed); 
                if (is.na(ii)) 0 else - sum(d * WintF.sl.alph * Xtime.deriv[, ii] - kii, na.rm = TRUE)},
            "both" = {ii <- match(i, indFixed); 
                - sum(d * (WintF.vl.alph * Xtime[, i] + 
                    WintF.sl.alph * if (is.na(ii)) 0 else Xtime.deriv[, ii]) - kii, na.rm = TRUE)}
        )
    }    
    c(sc1 + sc2)
}
