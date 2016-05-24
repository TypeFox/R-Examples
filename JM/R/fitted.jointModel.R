fitted.jointModel <-
function (object, process = c("Longitudinal", "Event"), 
        type = c("Marginal", "Subject", "EventTime", "Slope"), scale = c("survival", 
        "cumulative-Hazard", "log-cumulative-Hazard"), M = 200, ...) {
    if (!inherits(object, "jointModel"))
        stop("Use only with 'jointModel' objects.\n")
    process <- match.arg(process)
    type <- match.arg(type)
    scale <- match.arg(scale)
    method <- object$method
    if (process == "Longitudinal") {
        fitY <- c(object$x$X %*% object$coefficients$betas)
        names(fitY) <- names(object$y$y)
        if (type == "Marginal") fitY else if (type == "Subject") fitY + object$EB$Zb 
        else if (type == "EventTime") {
            fitY <- c(object$x$X %*% object$coefficients$betas) + object$EB$Zb
            fitYEvent <- if (method == "Cox-PH-GH") {
                c(object$x$Xtime2 %*% object$coefficients$betas + object$EB$Ztime2b)
            } else {
                c(object$x$Xtime %*% object$coefficients$betas + object$EB$Ztimeb)
            }
            id <- object$id
            idT <- object$x$idT
            fitYEvent <- unlist(mapply("c", split(fitY, id), split(fitYEvent, idT)), 
                use.names = FALSE)
            times <- object$times
            Time <- exp(object$y$logT)
            ind <- unlist(mapply("c", split(times, id), split(Time, idT)), 
                use.names = FALSE)
            fitYEvent <- fitYEvent[ind != 0]
            names(fitYEvent) <- seq_along(fitYEvent)
            fitYEvent
        } else if (type == "Slope") {
            derivForm <- object$derivForm
            indFixed <- derivForm$indFixed
            indRandom <- derivForm$indRandom
            TermsX.deriv <- object$termsYx.deriv
            TermsZ.deriv <- object$termsYz.deriv
            mfX.deriv <- model.frame(TermsX.deriv, data = object$data)
            mfZ.deriv <- model.frame(TermsZ.deriv, data = object$data)
            X.deriv <- model.matrix(derivForm$fixed, mfX.deriv)
            Z.deriv <- model.matrix(derivForm$random, mfZ.deriv)
            betas <- object$coefficients$betas
            b <- ranef(object)
            id <- object$id
            ff <- c(X.deriv %*% betas[indFixed] + 
                rowSums(Z.deriv * b[id, indRandom, drop = FALSE]))
            names(ff) <- names(object$y$y)
            Xtime.deriv <- object$x$Xtime.deriv
            Ztime.derivc <- object$x$Ztime.deriv
            ffEvent <- c(Xtime.deriv %*% betas[indFixed] + 
                rowSums(Ztime.deriv * b[, indRandom, drop = FALSE]))
            id <- object$id
            idT <- object$x$idT
            ffEvent <- unlist(mapply("c", split(ff, id), split(ffEvent, idT)), 
                use.names = FALSE)
            times <- object$times
            Time <- exp(object$y$logT)
            ind <- unlist(mapply("c", split(times, id), split(Time, idT)), 
                use.names = FALSE)
            ffEvent <- ffEvent[ind != 0]
            names(ffEvent) <- seq_along(ffEvent)
            ffEvent
        }
    } else {
        W1 <- object$x$W
        if (is.null(object$x[["Xtime"]])) {
            object$x[["Xtime"]] <- model.matrix(object$formYx, 
                                                model.frame(object$termsYx, data = object$data.id))
            object$x[["Ztime"]] <- model.matrix(object$formYz, 
                                                model.frame(object$termsYz, data = object$data.id))
            object$EB[["Ztimeb"]] <- rowSums(object$x[["Ztime"]] * ranef(object))
            
        }
        Y <- if (type == "Marginal") {
            D <- object$coefficients$D
            diag.D <- ncol(D) == 1 & (ncz <- nrow(D)) > 1
            b <- mvrnorm(M, rep(0, ncz), if (diag.D) diag(c(D)) else D)
            if (method == "Cox-PH-GH") {
                Zb <- object$x$Ztime2 %*% t(b)
                c(object$x$Xtime2 %*% object$coefficients$betas) + Zb
            } else {
                Zb <- object$x$Ztime %*% t(b)
                c(object$x$Xtime %*% object$coefficients$betas) + Zb
            }
        } else {
            if (method == "Cox-PH-GH") {
                c(object$x$Xtime2 %*% object$coefficients$betas + object$EB$Ztime2b)
            } else {
                c(object$x$Xtime %*% object$coefficients$betas + object$EB$Ztimeb)
            }
        }
        gammas <- object$coefficients$gammas
        alpha <- object$coefficients$alpha
        logT <- object$y$logT
        parameterization <- object$parameterization
        fitT <- if (method == "Cox-PH-GH") {
            indT <- object$indexes$indT
            lambda0 <- object$coefficients$lambda0[, "basehaz"]
            eta.tw <- if (!is.null(W1)) as.vector(W1 %*% gammas) else 0
            Ztime2b <- if (type == "Marginal") {
                object$x$Ztime2 %*% t(b)
            } else {
                rowSums(object$x$Ztime2 * object$EB$post.b[indT, ])
            }
            Y2 <- c(object$x$Xtime2 %*% object$coefficients$betas) + Ztime2b
            eta.s <- object$coefficients$alpha * Y2
            if (type == "Marginal") {
                S <- matrix(0, length(logT), ncol(eta.s))
                S[unique(indT), ] <- rowsum(lambda0[object$indexes$ind.L1] * exp(eta.s), 
                    indT, reorder = FALSE)
            } else {
                S <- numeric(length(logT))
                S[unique(indT)] <- tapply(lambda0[object$indexes$ind.L1] * exp(eta.s), indT, sum)
            }
            Haz <- exp(eta.tw) * S
            switch(scale,
                "survival" = exp(- Haz),
                "cumulative-Hazard" = Haz,
                "log-cumulative-Hazard" = log(Haz))
        } else if (method == "weibull-PH-GH") {
            WW <- if (is.null(W1)) as.matrix(rep(1, length(logT))) else cbind(1, W1)
            eta.tw <- as.vector(WW %*% gammas)
            sigma.t <- object$coefficients$sigma.t
            P <- object$x$P
            log.st <- log(object$x$st)
            wk <- rep(object$x$wk, length(logT))
            id.GK <- rep(seq_along(logT), each = object$control$GKk)
            if (parameterization %in% c("value", "both")) {
                Zsb <- if (type == "Marginal")
                    object$x$Zs %*% t(b)
                else
                    rowSums(object$x$Zs * object$EB$post.b[id.GK, ])
                Ys <- c(object$x$Xs %*% object$coefficients$betas) + Zsb
                eta.s <- c(object$x$Ws.intF.vl %*% object$coefficients$alpha) * Ys
            }
            if (parameterization %in% c("slope", "both")) {
                Zsb <- if (type == "Marginal")
                    object$x$Zs.deriv %*% t(b[, object$derivForm$indRandom, drop = FALSE])
                else
                    rowSums(object$x$Zs.deriv * object$EB$post.b[id.GK, object$derivForm$indRandom, drop = FALSE])
                Ys.deriv <- c(object$x$Xs.deriv %*% 
                    object$coefficients$betas[object$derivForm$indFixed]) + Zsb
                eta.s <- if (parameterization == "both")
                    eta.s + c(object$x$Ws.intF.sl %*% object$coefficients$Dalpha) * Ys.deriv 
                else
                    c(object$x$Ws.intF.sl %*% object$coefficients$Dalpha) * Ys.deriv
            }
            Haz <- exp(eta.tw) * P * rowsum(wk * exp(log(sigma.t) + 
                (sigma.t - 1) * log.st + eta.s), id.GK, reorder = FALSE)
            switch(scale,
                "survival" = exp(- Haz),
                "cumulative-Hazard" = Haz,
                "log-cumulative-Hazard" = log(Haz))
        } else if (method == "weibull-AFT-GH") {
            WW <- if (is.null(W1)) as.matrix(rep(1, length(logT))) else cbind(1, W1)
            eta.tw <- as.vector(WW %*% gammas)
            sigma.t <- object$coefficients$sigma.t
            P <- object$x$P
            log.st <- log(object$x$st)
            wk <- rep(object$x$wk, length(logT))
            id.GK <- rep(seq_along(logT), each = object$control$GKk)
            if (parameterization %in% c("value", "both")) {
                Zsb <- if (type == "Marginal")
                    object$x$Zs %*% t(b)
                else
                    rowSums(object$x$Zs * object$EB$post.b[id.GK, ])
                Ys <- c(object$x$Xs %*% object$coefficients$betas) + Zsb
                eta.s <- c(object$x$Ws.intF.vl %*% object$coefficients$alpha) * Ys
            }
            if (parameterization %in% c("slope", "both")) {
                Zsb <- if (type == "Marginal")
                    object$x$Zs.deriv %*% t(b[, object$derivForm$indRandom, drop = FALSE])
                else
                    rowSums(object$x$Zs.deriv * object$EB$post.b[id.GK, object$derivForm$indRandom, drop = FALSE])
                Ys.deriv <- c(object$x$Xs.deriv %*% 
                    object$coefficients$betas[object$derivForm$indFixed]) + Zsb
                eta.s <- if (parameterization == "both")
                    eta.s + c(object$x$Ws.intF.sl %*% object$coefficients$Dalpha) * Ys.deriv 
                else
                    c(object$x$Ws.intF.sl %*% object$coefficients$Dalpha) * Ys.deriv
            }            
            Vi <- exp(eta.tw) * P * rowsum(wk * exp(eta.s), id.GK, reorder = FALSE); dimnames(Vi) <- NULL            
            Haz <- Vi^sigma.t
            switch(scale,
                "survival" = exp(- Haz),
                "cumulative-Hazard" = Haz,
                "log-cumulative-Hazard" = log(Haz))
        } else if (method == "spline-PH-GH") {
            eta.tw <- if (!is.null(W1)) as.vector(W1 %*% gammas) else 0
            P <- object$x$P
            wk <- rep(object$x$wk, length(logT))
            id.GK <- rep(seq_along(logT), each = object$control$GKk)
            if (parameterization %in% c("value", "both")) {
                Zsb <- if (type == "Marginal") {
                    object$x$Zs %*% t(b)
                } else {
                    bb <- object$EB$post.b[object$x$idT, , drop = FALSE]
                    rowSums(object$x$Zs * bb[id.GK, ])
                }
                Ys <- c(object$x$Xs %*% object$coefficients$betas) + Zsb
                eta.s <- c(object$x$Ws.intF.vl %*% object$coefficients$alpha) * Ys
            }
            if (parameterization %in% c("slope", "both")) {
                Zsb <- if (type == "Marginal") {
                    object$x$Zs.deriv %*% t(b[, object$derivForm$indRandom, drop = FALSE])
                } else {
                    bb <- object$EB$post.b[object$x$idT, , drop = FALSE]
                    rowSums(object$x$Zs.deriv * 
                        bb[id.GK, object$derivForm$indRandom, drop = FALSE])
                }
                Ys.deriv <- c(object$x$Xs.deriv %*% 
                    object$coefficients$betas[object$derivForm$indFixed]) + Zsb
                eta.s <- if (parameterization == "both")
                    eta.s + c(object$x$Ws.intF.sl %*% object$coefficients$Dalpha) * Ys.deriv 
                else
                    c(object$x$Ws.intF.sl %*% object$coefficients$Dalpha) * Ys.deriv
            }
            Haz <- exp(eta.tw) * P * rowsum(wk * exp(c(object$x$W2s %*% object$coefficients$gammas.bs) + eta.s), 
                id.GK, reorder = FALSE)
            switch(scale,
                "survival" = exp(- Haz),
                "cumulative-Hazard" = Haz,
                "log-cumulative-Hazard" = log(Haz))
        } else if (method == "piecewise-PH-GH") {
            WW <- W1
            eta.tw <- if (!is.null(WW)) as.vector(WW %*% gammas) else 0
            xi <- object$coefficients$xi
            nk <- object$control$GKk
            id.GK <- object$x$id.GK
            ind.K <- rep(unlist(lapply(object$y$ind.D, seq_len)), each = nk)
            wk <- unlist(lapply(object$y$ind.D, function (n) rep(object$x$wk, n)))
            P <- c(t(object$x$P))
            wkP <- wk * rep(P[!is.na(P)], each = nk)
            if (parameterization %in% c("value", "both")) {
                Zsb <- if (type == "Marginal")
                    object$x$Zs %*% t(b)
                else
                    rowSums(object$x$Zs * object$EB$post.b[id.GK, ])
                Ys <- c(object$x$Xs %*% object$coefficients$betas) + Zsb
                eta.s <- c(object$x$Ws.intF.vl %*% object$coefficients$alpha) * Ys
            }
            if (parameterization %in% c("slope", "both")) {
                Zsb <- if (type == "Marginal")
                    object$x$Zs.deriv %*% t(b[, object$derivForm$indRandom, drop = FALSE])
                else
                    rowSums(object$x$Zs.deriv * object$EB$post.b[id.GK, object$derivForm$indRandom, drop = FALSE])
                Ys.deriv <- c(object$x$Xs.deriv %*% 
                    object$coefficients$betas[object$derivForm$indFixed]) + Zsb
                eta.s <- if (parameterization == "both")
                    eta.s + c(object$x$Ws.intF.sl %*% object$coefficients$Dalpha) * Ys.deriv 
                else
                    c(object$x$Ws.intF.sl %*% object$coefficients$Dalpha) * Ys.deriv
            }
            Haz <- exp(eta.tw) * rowsum(xi[ind.K] * wkP * exp(eta.s), id.GK, reorder = FALSE)
            switch(scale,
                "survival" = exp(- Haz),
                "cumulative-Hazard" = Haz,
                "log-cumulative-Hazard" = log(Haz))
        } else {
            W2 <- splineDesign(object$knots, logT, ord = object$control$ord)
            WW <- if (is.null(W1)) W2 else cbind(W2, W1)
            eta <- c(WW %*% gammas) + Y * alpha
             switch(scale,
                "survival" = exp(- exp(eta)),
                "cumulative-Hazard" = exp(eta),
                "log-cumulative-Hazard" = eta)
        }
        fitT <- if (type == "Marginal") rowMeans(fitT) else c(fitT)
        names(fitT) <- names(logT)
        fitT
    }
}
