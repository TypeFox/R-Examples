residuals.jointModel <-
function (object, process = c("Longitudinal", "Event"), 
    type = c("Marginal", "Subject", "stand-Marginal", "stand-Subject", 
        "Martingale", "nullMartingale", "CoxSnell", "AFT"), MI = FALSE, 
    M = 50, time.points = NULL, return.data = FALSE, ...) {
    if (!inherits(object, "jointModel"))
        stop("Use only with 'jointModel' objects.\n")
    process <- match.arg(process)
    type <- match.arg(type)
    if (process == "Longitudinal") {
        # Observed data
        y <- object$y$y
        X <- object$x$X
        Z <- object$x$Z
        ncx <- ncol(X)
        ncz <- ncol(Z)
        id <- object$id
        # Fitted Values
        fitted.vals <- if (type == "Marginal" || type == "stand-Marginal") {
            as.vector(X %*% object$coefficients$betas)
        } else {
            as.vector(X %*% object$coefficients$betas + object$EB$Zb)
        }
        # Residuals
        resid.vals <- if (type == "Marginal" || type == "Subject") {
            as.vector(y - fitted.vals)
        } else if (type == "stand-Subject"){
            as.vector(y - fitted.vals) / object$coefficients$sigma
        } else {
            D <- object$coefficients$D
            if (nrow(D) != ncol(D))
                D <- diag(c(D))
            unlist(lapply(split(cbind(Z, as.vector(y - fitted.vals)), id), function (x) {
                M <- matrix(x, ncol = ncz + 1)
                z <- M[, - (ncz + 1), drop = FALSE]
                res <- M[, ncz + 1]
                out <- z %*% D %*% t(z)
                diag(out) <- diag(out) + object$coefficients$sigma^2
                solve(chol(out)) %*% res
            }))
        }
        if (!MI) {
            names(resid.vals) <- names(y)
            resid.vals
        } else {
            if (object$CompRisk)
                stop("residuals() with 'MI = TRUE' is not currently implemented for ",
                    "competing risks joint models.\n")
            logT <- object$y$logT
            d <- object$y$d
            Xtime <- object$x$Xtime
            Ztime <- object$x$Ztime
            Xtime.deriv <- object$x$Xtime.deriv
            Ztime.deriv <- object$x$Ztime.deriv            
            Xs <- object$x$Xs
            Zs <- object$x$Zs
            Xs.deriv <- object$x$Xs.deriv
            Zs.deriv <- object$x$Zs.deriv
            WintF.vl <- object$x$WintF.vl
            WintF.sl <- object$x$WintF.sl
            Ws.intF.vl <- object$x$Ws.intF.vl
            Ws.intF.sl <- object$x$Ws.intF.sl
            P <- object$x$P
            wk <- object$x$wk
            method <- object$method
            parameterization <- object$parameterization
            LongFormat <- object$LongFormat
            idT <- object$x$idT
            W1 <- object$x$W
            WW <- if (method == "Cox-PH-GH") {
                stop("multiple-imputation-based residuals are not available ", 
                    "for joint models with method = 'Cox-PH-GH'.\n")
            } else if (method == "piecewise-PH-GH") {
                ind.D <- object$y$ind.D
                nk <- object$control$GKk
                st <- object$x$st
                ind.K <- rep(unlist(lapply(ind.D, seq_len)), each = nk)
                wk <- unlist(lapply(ind.D, function (n) rep(object$x$wk, n)))
                wkP <- wk * rep(object$x$P, each = nk)
                W1
            } else if (method == "weibull-PH-GH" || method == "weibull-AFT-GH") {
                log.st <- log(object$x$st)
                if (is.null(W1)) as.matrix(rep(1, length(logT))) else cbind(1, W1)
            } else if (method == "spline-PH-GH") {
                W2 <- object$x$W2
                W2s <- object$x$W2s
                W1
            } else {
                W2 <- splineDesign(object$knots, logT, ord = object$control$ord)
                nk <- ncol(W2) 
                if (is.null(W1)) W2 else cbind(W2, W1)
            }
            ncww <- if (is.null(WW)) 0 else ncol(WW)
            n <- length(logT)
            ni <- as.vector(tapply(id, id, length))
            obs.times <- if (!object$timeVar %in% names(object$data)) {
                if (is.null(ot <- attr(time.points, "obs.times")))
                    stop("could not extract observed times from either the design ", 
                        "matrix for the longitudinal measurements or\n\tthe 'time.points' argument.\n")
                else
                    ot
            } else { 
                object$data[[object$timeVar]]
            }
            environment(MI.fixed.times) <- environment(MI.random.times) <- environment()
            if (inherits(time.points, "weibull.frailty")) {
                MI.random.times(time.points)
            } else {
                MI.fixed.times(time.points)
            }
        }
    } else {
        #fits <- fitted(object, process = "Event", type = "Subject", scale = "cumulative-Hazard")
        #events <- object$y$d
        fits <- cumHaz(object, type == "nullMartingale")
        ni <- tapply(object$id, object$id, length)
        events <- rep(object$y$d, ni)
        events <- ave(events, object$id, FUN = function (x) c(rep(0, length(x)-1), x[1]))
        if (type == "AFT") {
            if (object$method == "weibull-AFT-GH") {
                fitted(object, process = "Event", type = "Subject", scale = "log-cumulative-Hazard")
            } else {
                warning("AFT residuals are only calculated for the Weibull AFT model; ",
                    "martingale residuals are calculated instead.\n")
                events - fits
            }
        } else if (type == "CoxSnell") {
            if (object$method %in% c("weibull-PH-GH", "weibull-AFT-GH", 
                "piecewise-PH-GH", "spline-PH-GH", "ch-Laplace")) {
                fitted(object, process = "Event", type = "Subject", scale = "cumulative-Hazard")
            } else {
                warning("CoxSnell residuals are only calculated for the parametric survival ",
                    "models; martingale residuals are calculated instead.\n")
                events - fits                
            }
        } else {
            events - fits
        }
    }
}
