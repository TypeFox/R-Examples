log.posterior.b2 <-
function (object) {
    method <- object$method
    timeVar <- object$timeVar
    interFact <- object$interFact
    parameterization <- object$parameterization
    derivForm <- object$derivForm
    indFixed <- derivForm$indFixed
    indRandom <- derivForm$indRandom
    n <- object$n
    id <- object$id
    idT <- object$x$idT
    y <- object$y$y
    ind.D <- object$y$ind.D
    logT <- object$y$logT
    d <- object$y$d
    strt <- object$y$strata
    X <- object$x$X
    Xs <- object$x$Xs
    Xs.deriv <- object$x$Xs.deriv
    Xtime <- object$x$Xtime
    Xtime.deriv <- object$x$Xtime.deriv
    Z <- object$x$Z
    Zs <- object$x$Zs
    Zs.deriv <- object$x$Zs.deriv
    Ztime <- object$x$Ztime
    Ztime.deriv <- object$x$Ztime.deriv
    W <- object$x$W
    W2 <- object$x$W2
    W2s <- object$x$W2s
    WintF.vl <- object$x$WintF.vl
    Ws.intF.vl <- object$x$Ws.intF.vl
    WintF.sl <- object$x$WintF.sl
    Ws.intF.sl <- object$x$Ws.intF.sl
    ####
    P <- object$x$P
    st <- object$x$st
    wk <- object$x$wk
    id.GK <- object$x$id.GK
    if (is.null(id.GK))
        id.GK <- rep(seq_along(logT), each = object$control$GKk)
    Q <- object$x$Q
    ncx <- ncol(X)
    ncz <- ncol(Z)
    ncww <- ncol(W)
    lag <- object$y$lag
    betas <- object$coefficients$betas
    sigma <- object$coefficients$sigma
    D <- object$coefficients$D
    diag.D <- ncol(D) == 1 & nrow(D) > 1
    D <- if (diag.D) diag(c(D)) else D
    gammas <- object$coefficients$gammas
    alpha <- object$coefficients$alpha
    Dalpha <- object$coefficients$Dalpha
    sigma.t <- object$coefficients$sigma.t
    xi <- object$coefficients$xi
    gammas.bs <- object$coefficients$gammas.bs
    ####
    ff <- function (b, i) {
        id.i <- id %in% i
        idT.i <- idT %in% i
        id.GK.i <- id.GK %in% which(idT.i)
        X.i <- X[id.i, , drop = FALSE]
        Z.i <- Z[id.i, , drop = FALSE]
        mu.y <- as.vector(X.i %*% betas) + rowSums(Z.i * rep(b, each = nrow(Z.i)))
        logNorm <- dnorm(y[id.i], mu.y, sigma, TRUE)
        log.p.yb <- sum(logNorm)
        log.p.b <- dmvnorm(b, rep(0, ncz), D, TRUE)
        eta.tw <- if (!is.null(W)) {
            if (method %in% c("weibull-PH-GH", "weibull-AFT-GH"))
                as.vector(cbind(1, W[idT.i, , drop = FALSE]) %*% gammas)
            else
                as.vector(W[idT.i, , drop = FALSE] %*% gammas)
        } else {
            0
        }
        if (parameterization %in% c("value", "both")) {
            Xtime.i <- Xtime[idT.i, , drop = FALSE]
            Ztime.i <- Ztime[idT.i, , drop = FALSE]
            Y <- c(Xtime.i %*% betas) + rowSums(Ztime.i * rep(b, each = nrow(Ztime.i)))
            Xs.i <- Xs[id.GK.i, , drop = FALSE]
            Zs.i <- Zs[id.GK.i, , drop = FALSE]
            Ys <- as.vector(Xs.i %*% betas + rowSums(Zs.i * rep(b, each = nrow(Zs.i))))
            WintF.vl.i <- WintF.vl[idT.i, , drop = FALSE]
            Ws.intF.vl.i <- Ws.intF.vl[id.GK.i, , drop = FALSE]
        }
        if (parameterization %in% c("slope", "both")) {
            Xtime.deriv.i <- Xtime.deriv[idT.i, , drop = FALSE]
            Ztime.deriv.i <- Ztime.deriv[idT.i, , drop = FALSE]
            Y.deriv <- c(Xtime.deriv.i %*% betas[indFixed]) + 
                rowSums(Ztime.deriv.i * rep(b[indRandom], each = nrow(Ztime.deriv.i)))
            Xs.deriv.i <- Xs.deriv[id.GK.i, , drop = FALSE]
            Zs.deriv.i <- Zs.deriv[id.GK.i, , drop = FALSE]
            Ys.deriv <- as.vector(Xs.deriv.i %*% betas[indFixed]) + 
                rowSums(Zs.deriv.i * rep(b[indRandom], each = nrow(Zs.deriv.i)))
            WintF.sl.i <- WintF.sl[idT.i, , drop = FALSE]
            Ws.intF.sl.i <- Ws.intF.sl[id.GK.i, , drop = FALSE]
        }
        tt <- switch(parameterization,
            "value" = c(WintF.vl.i %*% alpha) * Y, 
            "slope" = c(WintF.sl.i %*% Dalpha) * Y.deriv,
            "both" = c(WintF.vl.i %*% alpha) * Y + 
                c(WintF.sl.i %*% Dalpha) * Y.deriv)
        ss <- switch(parameterization,
            "value" = c(Ws.intF.vl.i %*% alpha) * Ys, 
            "slope" = c(Ws.intF.sl.i %*% Dalpha) * Ys.deriv,
            "both" = c(Ws.intF.vl.i %*% alpha) * Ys + 
                c(Ws.intF.sl.i %*% Dalpha) * Ys.deriv)
        ###
        log.p.tb <- if (method == "weibull-PH-GH") {
            st.i <- st[id.GK.i]
            Vi <- exp(log(sigma.t) + (sigma.t - 1) * log(st.i) + ss)
            log.S <- - exp(eta.tw) * P[i] * sum(wk * Vi)
            log.h <- log(sigma.t) + (sigma.t - 1) * logT[i] + eta.tw + tt
            d[i] * log.h + log.S
        } else if (method == "weibull-AFT-GH") {
            Vi <- exp(eta.tw) * P[i] * sum(wk * exp(ss))
            log.S <- - Vi^sigma.t
            log.h <- log(sigma.t) + (sigma.t - 1) * log(Vi) + eta.tw + tt
            d[i] * log.h + log.S
        } else if (method == "spline-PH-GH") {
            W2.i <- W2[idT.i, , drop = FALSE]
            W2s.i <- W2s[id.GK.i, , drop = FALSE]
            Vi <- exp(c(W2s.i %*% gammas.bs) + ss)
            kk.i <- sum(idT.i)
            log.S <- - exp(eta.tw) * P[idT.i] * 
                tapply(rep(wk, kk.i) * Vi, 
                    rep(seq_len(kk.i), each = length(wk)), sum)
            log.h <- eta.tw + tt + c(W2.i %*% gammas.bs) 
            sum(d[idT.i] * log.h + log.S)
        } else if (method == "piecewise-PH-GH") {
            nk <- object$control$GKk
            wkP <- rep(wk, sum(id.GK.i)/nk) * rep(P, each = nk)[id.GK.i]
            log.S <- - exp(eta.tw) * sum(xi[ind.D[i]] * wkP * exp(ss))
            log.h <- log(xi[ind.D[i]]) + eta.tw + tt
            as.vector(d[i] * log.h + log.S)
        }
        -(log.p.yb + log.p.tb + log.p.b)
    }
    modes <- ranef(object)
    hessians <- vector("list", n)
    for (i in seq_len(n)) {
        opt <- optim(modes[i, ], ff, i = i, method = "BFGS", hessian = TRUE)
        modes[i, ] <- opt$par
        hessians[[i]] <- solve(opt$hessian)
    }
    rownames(modes) <- if (!object$CompRisk && !object$LongFormat) {
        names(object$y$logT)
    } else {
        seq_len(nrow(modes))
    }
    names(hessians) <- rownames(modes)
    list(modes = modes, vars = hessians)
}
