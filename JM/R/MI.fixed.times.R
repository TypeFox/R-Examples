MI.fixed.times <-
function (time.points) {
    indFixed <- object$derivForm$indFixed
    indRandom <- object$derivForm$indRandom
    # indexes for missing data
    if (is.null(time.points))
        time.points <- sort(unique(obs.times))
    ind.miss <- id %in% which(ni < length(time.points))
    id.miss <- id[ind.miss]
    unq.id.miss <- unique(id.miss)
    ni.miss <- length(time.points) - ni
    ni.miss <- ni.miss[ni.miss > 0]
    id2.miss <- rep(seq_along(ni.miss), ni.miss)
    ni <- as.vector(tapply(id.miss, id.miss, length))
    id3.miss <- rep(seq_along(ni), ni)
    id.GK <- if (object$method %in% c("weibull-PH-GH", "weibull-AFT-GH", "spline-PH-GH")) {
        rep(idT %in% unq.id.miss, each = object$control$GKk)
    } else if (object$method == "piecewise-PH-GH") {
        object$x$id.GK %in% unq.id.miss
    } else NULL
    # observed data corresponding to patients with
    # one or more missing values
    y.missO <- y[ind.miss]
    logT.missO <- logT[unq.id.miss]
    d.missO <- d[unq.id.miss]
    X.missO <- X[ind.miss, ]
    Z.missO <- Z[ind.miss, ]
    keep <- idT %in% unq.id.miss
    idT.missO <- if (LongFormat) idT[keep] else idT[unq.id.miss]
    if (parameterization %in% c("value", "both")) {
        Xtime.missO <- Xtime[keep, , drop = FALSE]
        Ztime.missO <- Ztime[keep, , drop = FALSE]
        WintF.vl.missO <- WintF.vl[keep, , drop = FALSE]
    }
    if (parameterization %in% c("slope", "both")) {
        Xtime.deriv.missO <- Xtime.deriv[keep, , drop = FALSE]
        Ztime.deriv.missO <- Ztime.deriv[keep, , drop = FALSE]
        WintF.sl.missO <- WintF.sl[keep, , drop = FALSE]
    }
    if (method %in% c("weibull-PH-GH", "weibull-AFT-GH")) {
        P.missO <- P[unq.id.miss]
        log.st.missO <- log.st[id.GK]
        if (parameterization %in% c("value", "both")) {
            Xs.missO <- Xs[id.GK, , drop = FALSE]
            Zs.missO <- Zs[id.GK, , drop = FALSE]
        }
        if (parameterization %in% c("slope", "both")) {
            Xs.deriv.missO <- Xs.deriv[id.GK, , drop = FALSE]
            Zs.deriv.missO <- Zs.deriv[id.GK, , drop = FALSE]
        }
        Ws.intF.vl.missO <- Ws.intF.vl[id.GK, , drop = FALSE]
        Ws.intF.sl.missO <- Ws.intF.sl[id.GK, , drop = FALSE]
    }
    if (method == "spline-PH-GH") {
        P.missO <- P[keep]
        if (parameterization %in% c("value", "both")) {
            Xs.missO <- Xs[id.GK, , drop = FALSE]
            Zs.missO <- Zs[id.GK, , drop = FALSE]
        }
        if (parameterization %in% c("slope", "both")) {
            Xs.deriv.missO <- Xs.deriv[id.GK, , drop = FALSE]
            Zs.deriv.missO <- Zs.deriv[id.GK, , drop = FALSE]
        }
        W2s.missO <- W2s[id.GK, , drop = FALSE] 
        W2.missO <- W2[keep, , drop = FALSE]
        Ws.intF.vl.missO <- Ws.intF.vl[id.GK, , drop = FALSE]
        Ws.intF.sl.missO <- Ws.intF.sl[id.GK, , drop = FALSE]
    }
    if (method == "piecewise-PH-GH") {
        st.missO <- st[id.GK]
        ind.D.missO <- ind.D[unq.id.miss]
        ind.K.missO <- ind.K[id.GK]
        wkP.missO <- wkP[id.GK]
        if (parameterization %in% c("value", "both")) {
            Xs.missO <- Xs[id.GK, , drop = FALSE]
            Zs.missO <- Zs[id.GK, , drop = FALSE]
        }
        if (parameterization %in% c("slope", "both")) {
            Xs.deriv.missO <- Xs.deriv[id.GK, , drop = FALSE]
            Zs.deriv.missO <- Zs.deriv[id.GK, , drop = FALSE]
        }
        Ws.intF.vl.missO <- Ws.intF.vl[id.GK, , drop = FALSE]
        Ws.intF.sl.missO <- Ws.intF.sl[id.GK, , drop = FALSE]
    }
    WW.missO <- WW[unq.id.miss, , drop = FALSE]
    # observed data corresponding to patients with
    # one or more missing values
    mis.times <- unlist(lapply(split(obs.times, id), 
        function (x) time.points[!time.points %in% x]))
    dataM <- object$data[unq.id.miss, ]
    dataM <- dataM[id2.miss, ]
    dataM[[object$timeVar]] <- pmax(mis.times - object$y$lag, 0)
    mfX <- model.frame(object$termsYx, data = dataM)
    mfZ <- model.frame(object$termsYz, data = dataM)
    X.missM <- model.matrix(object$formYx, mfX)
    Z.missM <- model.matrix(object$formYz, mfZ)
    N.missM <- nrow(X.missM)
    n.missO <- length(unq.id.miss)
    # Estimated MLEs
    D <- object$coefficients$D
    diag.D <- ncz != ncol(D)
    list.thetas <- if (object$method == "weibull-PH-GH" || 
        object$method == "weibull-AFT-GH") {
        list(betas = object$coefficients$betas, 
            log.sigma = log(object$coefficients$sigma), 
            gammas = object$coefficients$gammas, 
            alpha = object$coefficients$alpha,
            Dalpha = object$coefficients$Dalpha, 
            log.sigma.t = log(object$coefficients$sigma.t),
            D = if (diag.D) log(D) else chol.transf(D))
    } else if (object$method == "spline-PH-GH") {
        list(betas = object$coefficients$betas, 
            log.sigma = log(object$coefficients$sigma),
            gammas = object$coefficients$gammas, 
            alpha = object$coefficients$alpha, 
            Dalpha = object$coefficients$Dalpha,
            gammas.bs = object$coefficients$gammas.bs, 
            D = if (diag.D) log(D) else chol.transf(D))
    } else if (object$method == "piecewise-PH-GH") {
        list(betas = object$coefficients$betas, 
            log.sigma = log(object$coefficients$sigma),
            gammas = object$coefficients$gammas,
            alpha = object$coefficients$alpha,
            Dalpha = object$coefficients$Dalpha,
            log.xi = log(object$coefficients$xi),
            D = if (diag.D) log(D) else chol.transf(D))
    }
    if (!is.null(object$scaleWB))
        list.thetas$log.sigma.t <- NULL
    list.thetas <- list.thetas[!sapply(list.thetas, is.null)]
    thetas <- unlist(as.relistable(list.thetas))
    V.thetas <- vcov(object)
    EBs <- ranef(object, postVar = TRUE)
    Var <- attr(EBs, "postVar")[unq.id.miss]
    EBs <- proposed.b <- EBs[unq.id.miss, , drop = FALSE]
    # Fitted values for corresponding to Y_i^m
    fitted.valsM <- if (type == "Marginal" || type == "stand-Marginal") {
        as.vector(X.missM %*% object$coefficients$betas)
    } else {
        as.vector(X.missM %*% object$coefficients$betas + 
            rowSums(Z.missM * EBs[id2.miss, , drop = FALSE]))
    }   
    current.b <- b.new <- EBs
    resid.valsM <- matrix(0, N.missM, M)
    environment(posterior.b) <- environment()
    old <- options(warn = (-1))
    on.exit(options(old))    
    for (m in 1:M) {
        # Step1: simulate new parameter values from a multivariate normal
        thetas.new <- mvrnorm(1, thetas, V.thetas)
        thetas.new <- relist(thetas.new, skeleton = list.thetas)
        betas.new <- thetas.new$betas
        sigma.new <- exp(thetas.new$log.sigma)
        gammas.new <- thetas.new$gammas
        alpha.new <- thetas.new$alpha
        Dalpha.new <- thetas.new$Dalpha
        D.new <- thetas.new$D
        D.new <- if (diag.D) exp(D.new) else chol.transf(D.new)
        if (object$method == "weibull-PH-GH" || object$method == "weibull-AFT-GH")
            sigma.t.new <- if (is.null(object$scaleWB)) 
                exp(thetas.new$log.sigma.t)
            else
                object$scaleWB
        if (object$method == "spline-PH-GH")
            gammas.bs.new <- thetas.new$gammas.bs
        if (object$method == "piecewise-PH-GH")
            xi.new <- exp(thetas.new$log.xi); Q <- object$x$Q
        # Step2: Simulate new values for the random effects
        eta.yx <- as.vector(X.missO %*% betas.new)
        #eta.yxT <- as.vector(Xtime.missO %*% betas.new)
        eta.tw <- if (!is.null(WW)) as.vector(WW.missO %*% gammas.new) else 0
        dmvt.current <- dmvt.proposed <- numeric(n.missO)
        for (i in 1:n.missO) {
            proposed.b[i, ] <- rmvt(1, EBs[i, ], Var[[i]], 4)
            tt <- dmvt(rbind(current.b[i, ], proposed.b[i, ]), 
                EBs[i, ], Var[[i]], 4, TRUE)
            dmvt.current[i] <- tt[1]
            dmvt.proposed[i] <- tt[2]
        }
        a <- pmin(exp(posterior.b(proposed.b) + dmvt.current - 
            posterior.b(current.b) - dmvt.proposed), 1)
        ind <- runif(n.missO) <= a
        b.new[ind, ] <- proposed.b[ind, ]
        current.b <- b.new
        # Step3: Simulate new Y_i^m and calculate residuals
        mu <- as.vector(X.missM %*% betas.new + 
            rowSums(Z.missM * b.new[id2.miss, , drop = FALSE]))
        y.new <- rnorm(N.missM, mu, sigma.new)
        resid.valsM[, m] <- y.new - fitted.valsM
    }
    mean.resid.valsM <- rowMeans(resid.valsM)
    if (type == "stand-Subject") {
        var.resid.valsM <- object$coefficients$sigma^2 + 
            apply(resid.valsM, 1, var)
        mean.resid.valsM <- mean.resid.valsM / sqrt(var.resid.valsM)
    }
    if (type == "stand-Marginal") {
        mean.resid.valsM <- unlist(lapply(split(cbind(Z.missM, 
            resid.valsM), id2.miss), function (x) {
            MM <- matrix(x, ncol = ncz + M)
            z <- MM[, 1:ncz, drop = FALSE]
            res <- MM[, -(1:ncz), drop = FALSE]
            V1 <- z %*% D %*% t(z)
            diag(V1) <- diag(V1) + object$coefficients$sigma^2
            rr <- res - rowMeans(res)
            V2 <- apply(rr, 2, function (y) y %o%y)
            V2 <- if (is.matrix(V2)) rowSums(V2) / (M - 1) else sum(V2) / (M - 1)
            dim(V2) <- c(nrow(rr), nrow(rr))
            solve(chol(V1 + V2)) %*% rowMeans(res)
            }))
        }
        resid.valsM <- apply(resid.valsM, 2, function (x) {
            if (type == "stand-Subject")
                x <- x / object$coefficients$sigma
            if (type == "stand-Marginal") {
                x <- unlist(lapply(split(cbind(Z.missM, x), id2.miss), function (y) {
                    M <- matrix(y, ncol = ncz + 1)
                    z <- M[, - (ncz + 1), drop = FALSE]
                    res <- M[, ncz + 1]
                    out <- z %*% D %*% t(z)
                    diag(out) <- diag(out) + object$coefficients$sigma^2
                    solve(chol(out)) %*% res
                }))
            }
            x
        })
        names(resid.vals) <- names(fitted.vals) <- names(y)
        names(fitted.valsM) <- names(mean.resid.valsM) <- rownames(resid.valsM) <- paste("m", 1:length(fitted.valsM), sep = "")
        list("fitted.values" = fitted.vals, "residuals" = resid.vals, "fitted.valsM" = fitted.valsM, 
             "mean.resid.valsM" = mean.resid.valsM, "resid.valsM" = resid.valsM, 
             "dataM" = if (return.data) dataM else NULL)
}
