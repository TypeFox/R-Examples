MI.random.times <-
function (time.points) {
    indFixed <- object$derivForm$indFixed
    indRandom <- object$derivForm$indRandom
    # indexes for missing data
    t.max <- if (is.null(tt <- attr(time.points, "t.max"))) max(obs.times) else tt
    max.visits <- if (is.null(tt <- attr(time.points, "max.visits"))) max(ni) * 5 else tt
    id.GK <- rep(TRUE, length(object$x$id.GK))
    y.missO <- y
    logT.missO <- logT
    d.missO <- d
    X.missO <- X
    Z.missO <- Z
    idT.missO <- object$x$idT
    if (parameterization %in% c("value", "both")) {
        Xtime.missO <- Xtime
        Ztime.missO <- Ztime
        WintF.vl.missO <- WintF.vl
        Ws.intF.vl.missO <- Ws.intF.vl
    }
    if (parameterization %in% c("slope", "both")) {
        Xtime.deriv.missO <- Xtime.deriv
        Ztime.deriv.missO <- Ztime.deriv
        WintF.sl.missO <- WintF.sl
        Ws.intF.sl.missO <- Ws.intF.sl
    }
    if (method %in% c("weibull-PH-GH", "weibull-AFT-GH")) {
        P.missO <- P
        log.st.missO <- log.st
        if (parameterization %in% c("value", "both")) {
            Xs.missO <- Xs
            Zs.missO <- Zs
        }
        if (parameterization %in% c("slope", "both")) {
            Xs.deriv.missO <- Xs.deriv
            Zs.deriv.missO <- Zs.deriv
        }
    }
    if (method == "spline-PH-GH") {
        P.missO <- P
        if (parameterization %in% c("value", "both")) {
            Xs.missO <- Xs
            Zs.missO <- Zs
        }
        if (parameterization %in% c("slope", "both")) {
            Xs.deriv.missO <- Xs.deriv
            Zs.deriv.missO <- Zs.deriv
        }
        W2s.missO <- W2s
        W2.missO <- W2
    }
    if (method == "piecewise-PH-GH") {
        st.missO <- st
        ind.D.missO <- ind.D
        ind.K.missO <- ind.K
        wkP.missO <- wkP
        if (parameterization %in% c("value", "both")) {
            Xs.missO <- Xs
            Zs.missO <- Zs
        }
        if (parameterization %in% c("slope", "both")) {
            Xs.deriv.missO <- Xs.deriv
            Zs.deriv.missO <- Zs.deriv
        }
    }
    WW.missO <- WW
    n.missO <- length(unique(idT.missO))
    id.miss <- id3.miss <- id
    # Estimated MLEs (Joint Model)
    D <- object$coefficients$D
    diag.D <- ncz != ncol(D)
    list.thetas <- if (object$method == "weibull-PH-GH" || object$method == "weibull-AFT-GH") {
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
    Var <- attr(EBs, "postVar")
    EBs <- proposed.b <- EBs
    # Visit Times
    U <- time.points$y[, 1]
    X.vs <- time.points$x
    ncx.vs <- ncol(X.vs)
    id.onevisit <- as.vector(which(tapply(id, id, length) == 1))
    id.mrvisits <- as.vector(which(tapply(id, id, length) > 1))
    ev.vs <- tapply(id, id, length) - 1
    ev.vs <- as.vector(ev.vs[ev.vs > 0])
    id.vs.fl <- rep(id.mrvisits, ev.vs)
    n.vs.one <- length(id.onevisit)
    n.vs.more <- length(id.mrvisits)
    n.vs <- n.vs.one + n.vs.more
    # Estimated MLEs (Visiting Model)
    thetas.vs <- c(time.points$coefficients$betas, 
        log(time.points$coefficients$scale), 
        log(time.points$coefficients$shape), 
        log(time.points$coefficients$var.frailty))
    Var.vs <- vcov(time.points)
    p.vs <- length(thetas.vs)
    # current value for random effects
    current.b <- b.new <- EBs
    environment(posterior.b) <- environment()
    fitted.valsM.lis <- resid.valsM.lis <- vector("list", M)
    old <- options(warn = (-1))
    on.exit(options(old))
    for (m in 1:M) {
        curr.y <- tapply(object$y$y, object$id, function (x) x[length(x)])
        new.visit <- last.visit <- tapply(obs.times, object$id, function (x) x[length(x)])
        # Step1: simulate new parameter values from a multivariate normal
        # joint model
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
            sigma.t.new <- if (is.null(object$scaleWB)) {
                exp(thetas.new$log.sigma.t)
            } else {
                object$scaleWB
            }
        if (object$method == "spline-PH-GH")
            gammas.bs.new <- thetas.new$gammas.bs
        if (object$method == "piecewise-PH-GH")
            xi.new <- exp(thetas.new$log.xi); Q <- object$x$Q
        # visiting model
        thetas.vs.new <- mvrnorm(1, thetas.vs, Var.vs)
        betas.vs.new <- thetas.vs.new[seq_len(ncx.vs)]
        scale.vs.new <- exp(thetas.vs.new[ncx.vs + 1])
        shape.vs.new <- exp(thetas.vs.new[ncx.vs + 2])
        var.fr.new <- exp(thetas.vs.new[ncx.vs + 3])
        # Step2: Simulate new values for the random effects
        # joint model
        eta.yx <- as.vector(X.missO %*% betas.new)
        eta.yxT <- as.vector(Xtime.missO %*% betas.new)
        eta.tw <- as.vector(WW.missO %*% gammas.new)
        dmvt.current <- dmvt.proposed <- numeric(n.missO)
        for (i in seq_len(n.missO)) {
            proposed.b[i, ] <- rmvt(1, EBs[i, ], Var[[i]], 4)
            tt <- dmvt(rbind(current.b[i, ], proposed.b[i, ]), EBs[i, ], Var[[i]], 4, TRUE)
            dmvt.current[i] <- tt[1]
            dmvt.proposed[i] <- tt[2]
        }
        a <- pmin(exp(posterior.b(proposed.b) + dmvt.current - posterior.b(current.b) - dmvt.proposed), 1)
        ind <- runif(n.missO) <= a
        b.new[ind, ] <- proposed.b[ind, ]
        current.b <- b.new
        # visiting model
        omega.new <- numeric(n)
        omega.new[id.onevisit] <- rgamma(n.vs.one, 1/var.fr.new, 1/var.fr.new)
        exp.eta.vs <- exp(X.vs %*% betas.vs.new)
        omega.new[id.mrvisits] <- rgamma(n.vs.more, 1/var.fr.new + ev.vs, 
            1/var.fr.new + scale.vs.new * tapply(c(U^shape.vs.new * exp.eta.vs), id.vs.fl, sum))
        # Step2: Simulate new values for responses
        fitted.valsM <- resid.valsM <- Visit.Times <- matrix(as.numeric(NA), n, max.visits)
        Z.missM.lis <- vector("list", max.visits)
        ii <- 1
        while (any(new.visit[!is.na(new.visit)] < t.max)) {
            data.vs <- time.points$data[!duplicated(time.points$data[[time.points$nam.id]]), ]
            if (!is.null(nam <- attr(time.points, "prev.y")))
                data.vs[[nam]] <- curr.y
            mf.vs <- model.frame(time.points$terms, data = data.vs, na.action = NULL)
            X.vs.new <- model.matrix(formula(time.points), mf.vs)[, -1, drop = FALSE]
            mu.vs <- c(log(scale.vs.new) + X.vs.new %*% betas.vs.new) + log(omega.new) / shape.vs.new
            u.new <- rweibull(n, shape.vs.new, 1 / exp(mu.vs))
            Visit.Times[, ii] <- new.visit <- last.visit + u.new
            ind.tmax <- new.visit > t.max
            dataM <- object$data.id
            dataM[object$timeVar] <- pmax(new.visit - object$y$lag, 0)
            mfX <- model.frame(object$termsYx, data = dataM, na.action = NULL)
            mfZ <- model.frame(object$termsYz, data = dataM, na.action = NULL)
            X.missM <- model.matrix(object$formYx, mfX)
            Z.missM.lis[[ii]] <- Z.missM <- model.matrix(object$formYz, mfZ)
            fitted.valsM[, ii] <- if (type == "Marginal" || type == "stand-Marginal") {
                as.vector(X.missM %*% object$coefficients$betas)
            } else {
                as.vector(X.missM %*% object$coefficients$betas + rowSums(Z.missM * b.new))
            }
            mu <- as.vector(X.missM %*% betas.new + rowSums(Z.missM * b.new))
            y.new <- rnorm(n, mu, sigma.new)
            resid.valsM[, ii] <- y.new - fitted.valsM[, ii]
            Visit.Times[ind.tmax, ii] <- fitted.valsM[ind.tmax, ii] <- resid.valsM[ind.tmax, ii] <- as.numeric(NA)
            curr.y <- y.new
            last.visit <- new.visit
            ii <- ii + 1
            if (ii > max.visits)
                break
        }
        na.ind <- colSums(is.na(fitted.valsM)) != n
        Visit.Times <- Visit.Times[, na.ind]
        fitted.valsM <- fitted.valsM[, na.ind]
        resid.valsM <- resid.valsM[, na.ind]
        Z.missM <- do.call(rbind, Z.missM.lis[na.ind])
        id2.miss <- rep(1:n, ncol(resid.valsM))
        if (type == "stand-Subject")
            resid.valsM <- resid.valsM / object$coefficients$sigma
        if (type == "stand-Marginal") {
            resid.valsM <- unlist(lapply(split(cbind(Z.missM, c(resid.valsM)), id2.miss), function (y) {
                M <- matrix(y, ncol = ncz + 1)
                z <- M[, - (ncz + 1), drop = FALSE]
                res <- M[, ncz + 1]
                nz <- nrow(M)
                result <- rep(as.numeric(NA), nz)
                na.ind <- !is.na(res)
                if (all(!na.ind)) {
                    result
                } else {
                    out <- z[na.ind, , drop = FALSE] %*% D %*% t(z[na.ind, , drop = FALSE])
                    diag(out) <- diag(out) + object$coefficients$sigma^2
                    result[na.ind] <- solve(chol(out)) %*% res[na.ind]
                    result
                }
            }))
        }
        fitted.valsM.lis[[m]] <- fitted.valsM
        resid.valsM.lis[[m]] <- resid.valsM
    }
    names(resid.vals) <- names(fitted.vals) <- names(y)
    names(fitted.valsM.lis) <- names(resid.valsM.lis) <- paste("m", seq_len(M), sep = "")
    fitted.valsM.lis <- lapply(fitted.valsM.lis, function (x) {
        dimnames(x) <- list(1:n, paste("time", seq_len(ncol(x)), sep = ""))
        x
    })
    resid.valsM.lis <- if (type == "stand-Marginal") {
        resid.valsM.lis
    } else {
        lapply(resid.valsM.lis, function (x) {
            dimnames(x) <- list(1:n, paste("time", seq_len(ncol(x)), sep = ""))
            x
        })
    }
    list("fitted.values" = fitted.vals, "residuals" = resid.vals, "fitted.valsM" = fitted.valsM.lis, 
         "mean.resid.valsM" = NULL, "resid.valsM" = resid.valsM.lis, 
         "dataM" = NULL)
}
