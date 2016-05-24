log.posterior.b <-
function (b, y, Mats, method, ii) {
    id.i <- id %in% ii
    idT.i <- idT %in% ii
    X.i <- X[id.i, , drop = FALSE]
    Z.i <- Z[id.i, , drop = FALSE]
    mu.y <- as.vector(X.i %*% betas.new) + rowSums(Z.i * rep(b, each = nrow(Z.i)))
    logNorm <- dnorm(y[id.i], mu.y, sigma.new, TRUE)
    log.p.yb <- sum(logNorm)
    log.p.b <- dmvnorm(b, rep(0, ncol(Z)), D.new, TRUE)
    st <- Mats[[ii]]$st
    wk <- Mats[[ii]]$wk
    P <- Mats[[ii]]$P
    Xs <- Mats[[ii]]$Xs
    Zs <- Mats[[ii]]$Zs
    Xs.deriv <- Mats[[ii]]$Xs.deriv
    Zs.deriv <- Mats[[ii]]$Zs.deriv
    Ws.intF.vl <- Mats[[ii]]$Ws.intF.vl
    Ws.intF.sl <- Mats[[ii]]$Ws.intF.sl
    ind <- Mats[[ii]]$ind
    if (parameterization %in% c("value", "both"))
        Ys <- as.vector(Xs %*% betas.new + rowSums(Zs * rep(b, each = nrow(Zs))))
    if (parameterization %in% c("slope", "both"))
        Ys.deriv <- as.vector(Xs.deriv %*% betas.new[indFixed]) + 
            rowSums(Zs.deriv * rep(b[indRandom], each = nrow(Zs)))
    tt <- switch(parameterization,
        "value" = c(Ws.intF.vl %*% alpha.new) * Ys, 
        "slope" = c(Ws.intF.sl %*% Dalpha.new) * Ys.deriv,
        "both" = c(Ws.intF.vl %*% alpha.new) * Ys + 
            c(Ws.intF.sl %*% Dalpha.new) * Ys.deriv)
    eta.tw <- if (!is.null(W)) {
        if (!LongFormat)
            as.vector(W[ii, , drop = FALSE] %*% gammas.new)
        else
            as.vector(W[idT.i, , drop = FALSE] %*% gammas.new)
    } else 0
    log.survival <- if (method == "weibull-PH-GH") {
        Vi <- exp(log(sigma.t.new) + (sigma.t.new - 1) * log(st) + tt)
        - exp(eta.tw) * P * sum(wk * Vi)
    } else if (method == "weibull-AFT-GH") {
        Vi <- exp(eta.tw) * P * sum(wk * exp(tt))
        - Vi^sigma.t.new
    } else if (method == "spline-PH-GH") {
        W2s <- if (length(kn <- object$control$knots) == 1) {
            splineDesign(unlist(kn, use.names = FALSE), st, 
                ord = object$control$ord, outer.ok = TRUE)
        } else {
            strt.i <- strt[ii]
            w2s <- lapply(kn, function (kn) 
                splineDesign(kn, st, ord = object$control$ord, outer.ok = TRUE))
            ll <- match(strt.i, names(w2s))
            w2s[-ll] <- lapply(w2s[-ll], function (m) {m[, ] <- 0; m})
            do.call(cbind, w2s)
        }
        Vi <- exp(c(W2s %*% gammas.bs.new) + tt)
        idT <- rep(seq_along(P), each = object$control$GKk)
        - sum(exp(eta.tw) * P * tapply(wk * Vi, idT, sum))
    } else if (method == "piecewise-PH-GH") {
        P <- P[!is.na(P)]
        ind.K <- rep(seq_len(ind), each = 7)
        wk <- rep(wk, ind)
        wkP <- wk * rep(P, each = 7)
        eta.tw <- if (!is.null(W)) as.vector(W[i, , drop = FALSE] %*% gammas.new) else 0 
        - exp(eta.tw) * sum(xi.new[ind.K] * wkP * exp(tt))
    }
    log.p.yb + log.survival + log.p.b
}
