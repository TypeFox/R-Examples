plot.jointModel <- function (x, which = 1:4, caption = c("Residuals vs Fitted", "Normal Q-Q", "Marginal Survival", 
    "Marginal Cumulative Hazard", "Marginal log Cumulative Hazard", "Baseline Hazard", "Cumulative Baseline Hazard", 
    "Subject-specific Survival", "Subject-specific Cumulative Hazard", "Subject-specific log Cumulative Hazard"), 
    survTimes = NULL, main = "", ask = prod(par("mfcol")) < length(which) && dev.interactive(), ..., ids = NULL, 
    add.smooth = getOption("add.smooth"), add.qqline = TRUE, add.KM = FALSE, cex.caption = 1, return = FALSE) {
    if (!inherits(x, "jointModel"))
        stop("Use only with 'jointModel' objects.\n")
    if (!is.numeric(which) || any(which < 1) || any(which > 10))
        stop("'which' must be in 1:10.\n")
    show <- rep(FALSE, 10)
    if ((x$LongFormat || x$CompRisk) && any(which %in% c(3:5, 8:10))) {
        warning("\n  the marginal and subject-specific survival and cumulative ", 
            "hazard functions\n  are not currently implemented for joint models ", 
            "with exogenous time-dependent\n  covariates or competing risks.")
        which <- 1:2
    }
    show[which] <- TRUE
    method <- x$method
    if (any(show[6], show[7]) && method != "Cox-PH-GH") {
        show[6] <- show[7] <- FALSE
        warning("the baseline hazard and the cumulative baseline hazard are only plotted for the 'Cox-PH-GH' method.\n")
    }   
    if (any(show[c(3:5, 8:10)])) {
        if (is.null(ids))
            ids <- seq_len(x$n)
        if (is.null(survTimes) || !is.numeric(survTimes))
            survTimes <- seq(min(exp(x$y$logT)), max(exp(x$y$logT)), length.out = 31)
        log.survTimes <- log(survTimes)
        nt <- length(survTimes)
        n <- x$n
        W1 <- x$x$W
        gammas <- x$coefficients$gammas
        alpha <- x$coefficients$alpha
        Dalpha <- x$coefficients$Dalpha
        parameterization <- x$parameterization
        derivForm <- x$derivForm
        indFixed <- derivForm$indFixed
        indRandom <- derivForm$indRandom
        fitT <- if (method == "Cox-PH-GH") {
            lambda0 <- x$coefficients$lambda0[, "basehaz"]
            unqT <- x$coefficients$lambda0[, "time"]
            T.mat <- matrix(exp(log.survTimes), nrow = n, ncol = nt, byrow = TRUE)
            eta.tw <- if (!is.null(W1)) as.vector(W1 %*% gammas) else rep(0, n)
            Haz <- matrix(0, n, nt)
            for (i in 1:nt) {
                times <- lapply(T.mat[, i], function (t) unqT[t >= unqT])
                ind.len <- sapply(times, length)
                indT <- rep(1:nrow(x$data.id), ind.len)
                data.id2 <- x$data.id[indT, ]
                data.id2[[x$timeVar]] <- pmax(unlist(times, use.names = FALSE) - x$y$lag, 0)
                mfX <- model.frame(x$termsYx, data = data.id2)
                mfZ <- model.frame(x$termsYz, data = data.id2)
                Xtime2 <- model.matrix(x$formYx, mfX)
                Ztime2 <- model.matrix(x$formYz, mfZ)
                nk <- as.vector(sapply(split(indT, indT), length))
                ind.L1 <- unlist(lapply(nk, seq, from = 1))
                Y2 <- c(Xtime2 %*% x$coefficients$betas + rowSums(Ztime2 * x$EB$post.b[indT, ]))
                eta.s <- alpha * Y2
                S <- numeric(n)
                S[unique(indT)] <- tapply(lambda0[ind.L1] * exp(eta.s), indT, sum)
                Haz[, i] <- exp(eta.tw) * S
            }
            list("survival" = exp(- Haz), "cumulative-Hazard" = Haz, "log-cumulative-Hazard" = log(Haz))
        } else if (method == "weibull-PH-GH") {
            T.mat <- matrix(exp(log.survTimes), nrow = n, ncol = nt, byrow = TRUE)
            WW <- if (is.null(W1)) as.matrix(rep(1, n)) else cbind(1, W1)
            eta.tw <- as.vector(WW %*% gammas)
            sigma.t <- x$coefficients$sigma.t
            b <- x$EB$post.b
            wk <- gaussKronrod(x$control$GKk)$wk
            sk <- gaussKronrod(x$control$GKk)$sk
            id.GK <- rep(seq_len(n), each = x$control$GKk)
            Haz <- matrix(0, n, nt)
            for (i in 1:nt) {
                P <- T.mat[, i] / 2
                st <- outer(P, sk + 1)
                data.id <- x$data.id[id.GK, ]
                data.id[[x$timeVar]] <- pmax(c(t(st)) - x$y$lag, 0)
                if (parameterization %in% c("value", "both")) {
                    mfX <- model.frame(x$termsYx, data = data.id)
                    mfZ <- model.frame(x$termsYz, data = data.id)
                    Xs <- model.matrix(x$formYx, mfX)
                    Zs <- model.matrix(x$formYz, mfZ)
                    Ys <- c(Xs %*% x$coefficients$betas) + rowSums(Zs * b[id.GK, , drop = FALSE])
                    eta.s <- c(x$x$WintF.vl[id.GK, , drop = FALSE] %*% alpha) * Ys
                }
                if (parameterization %in% c("slope", "both")) {
                    mfX.deriv <- model.frame(x$termsYx.deriv, data = data.id)
                    mfZ.deriv <- model.frame(x$termsYz.deriv, data = data.id)
                    Xs.deriv <- model.matrix(derivForm$fixed, mfX.deriv)
                    Zs.deriv <- model.matrix(derivForm$random, mfZ.deriv)
                    Ys.deriv <- c(Xs.deriv %*% x$coefficients$betas[indFixed]) +
                        if (length(indRandom) > 1 || indRandom) 
                            rowSums(Zs.deriv * b[id.GK, indRandom, drop = FALSE])
                        else
                            0
                    eta.s <- if (parameterization == "both")
                        eta.s + c(x$x$WintF.sl[id.GK, , drop = FALSE] %*% Dalpha) * Ys.deriv
                    else
                        c(x$x$WintF.sl[id.GK, , drop = FALSE] %*% Dalpha) * Ys.deriv
                }
                log.st <- log(c(t(st)))
                Haz[, i] <- exp(eta.tw) * P * rowsum(wk * exp(log(sigma.t) + (sigma.t - 1) * log.st + eta.s), 
                    id.GK, reorder = FALSE)
            }
            list("survival" = exp(- Haz),
                 "cumulative-Hazard" = Haz,
                 "log-cumulative-Hazard" = log(Haz))
        } else if (method == "weibull-AFT-GH") {
            T.mat <- matrix(exp(log.survTimes), nrow = n, ncol = nt, byrow = TRUE)
            WW <- if (is.null(W1)) as.matrix(rep(1, n)) else cbind(1, W1)
            eta.tw <- as.vector(WW %*% gammas)
            sigma.t <- x$coefficients$sigma.t
            b <- x$EB$post.b
            wk <- gaussKronrod(x$control$GKk)$wk
            sk <- gaussKronrod(x$control$GKk)$sk
            id.GK <- rep(seq_len(n), each = x$control$GKk)
            Haz <- matrix(0, n, nt)
            for (i in 1:nt) {
                P <- T.mat[, i] / 2
                st <- outer(P, sk + 1)
                data.id <- x$data.id[id.GK, ]
                data.id[[x$timeVar]] <- pmax(c(t(st)) - x$y$lag, 0)
                if (parameterization %in% c("value", "both")) {
                    mfX <- model.frame(x$termsYx, data = data.id)
                    mfZ <- model.frame(x$termsYz, data = data.id)
                    Xs <- model.matrix(x$formYx, mfX)
                    Zs <- model.matrix(x$formYz, mfZ)
                    Ys <- c(Xs %*% x$coefficients$betas) + rowSums(Zs * b[id.GK, , drop = FALSE])
                    eta.s <- c(x$x$WintF.vl[id.GK, , drop = FALSE] %*% alpha) * Ys
                }
                if (parameterization %in% c("slope", "both")) {
                    mfX.deriv <- model.frame(x$termsYx.deriv, data = data.id)
                    mfZ.deriv <- model.frame(x$termsYz.deriv, data = data.id)
                    Xs.deriv <- model.matrix(derivForm$fixed, mfX.deriv)
                    Zs.deriv <- model.matrix(derivForm$random, mfZ.deriv)
                    Ys.deriv <- c(Xs.deriv %*% x$coefficients$betas[indFixed]) +
                        if (length(indRandom) > 1 || indRandom) 
                            rowSums(Zs.deriv * b[id.GK, indRandom, drop = FALSE])
                        else
                            0
                    eta.s <- if (parameterization == "both")
                        eta.s + c(x$x$WintF.sl[id.GK, , drop = FALSE] %*% Dalpha) * Ys.deriv
                    else
                        c(x$x$WintF.sl[id.GK, , drop = FALSE] %*% Dalpha) * Ys.deriv
                }
                log.st <- log(c(t(st)))
                Vi <- exp(eta.tw) * P * rowsum(wk * exp(eta.s), id.GK, reorder = FALSE); dimnames(Vi) <- NULL
                Haz[, i] <- Vi^sigma.t
            }
            list("survival" = exp(- Haz),
                 "cumulative-Hazard" = Haz,
                 "log-cumulative-Hazard" = log(Haz))
        } else if (method == "piecewise-PH-GH") {
            T.mat <- matrix(exp(log.survTimes), nrow = n, ncol = nt, byrow = TRUE)
            Q <- x$x$Q
            qs <- c(0, x$control$knots, max(exp(x$y$logT)) + 1)
            WW <- W1
            eta.tw <- if (!is.null(WW)) as.vector(WW %*% gammas) else 0
            xi <- x$coefficients$xi
            b <- x$EB$post.b
            sk <- gaussKronrod(x$control$GKk)$sk
            nk <- length(sk)
            Haz <- matrix(0, n, nt)
            for (i in 1:nt) {
                ind.D <- findInterval(T.mat[, i], qs, rightmost.closed = TRUE)
                Tiq <- outer(T.mat[, i], qs, pmin)
                Lo <- Tiq[, 1:Q]
                Up <- Tiq[, 2:(Q+1)]
                T <- Up - Lo
                P <- T / 2
                P[P < x$control$tol3] <- as.numeric(NA)
                P1 <- (Up + Lo) / 2
                st <- matrix(0, n, nk*Q)
                skQ <- rep(sk, Q)
                for (ii in 1:n) {
                    st[ii, ] <- rep(P[ii, ], each = nk) * skQ + rep(P1[ii, ], each = nk)
                }
                data.id2 <- x$data.id[rep(1:n, each = nk*Q), ]
                data.id2[[x$timeVar]] <- pmax(c(t(st)) - x$y$lag, 0)
                id.GK <- rep(1:n, rowSums(!is.na(st)))
                if (parameterization %in% c("value", "both")) {
                    mfX <- model.frame(x$termsYx, data = data.id2)
                    mfZ <- model.frame(x$termsYz, data = data.id2)
                    Xs <- model.matrix(x$formYx, mfX)
                    Zs <- model.matrix(x$formYz, mfZ)
                    #Zs <- Zs[!is.na(data.id2[[x$timeVar]]), ]
                    Ys <- c(Xs %*% x$coefficients$betas) + rowSums(Zs * b[id.GK, , drop = FALSE])
                    eta.s <- c(x$x$WintF.vl[id.GK, , drop = FALSE] %*% alpha) * Ys
                }
                if (parameterization %in% c("slope", "both")) {
                    mfX.deriv <- model.frame(x$termsYx.deriv, data = data.id2)
                    mfZ.deriv <- model.frame(x$termsYz.deriv, data = data.id2)
                    Xs.deriv <- model.matrix(derivForm$fixed, mfX.deriv)
                    Zs.deriv <- model.matrix(derivForm$random, mfZ.deriv)
                    #Zs.deriv <- Zs.deriv[!is.na(data.id2[[x$timeVar]]), ]
                    Ys.deriv <- c(Xs.deriv %*% x$coefficients$betas[indFixed]) +
                        if (length(indRandom) > 1 || indRandom) 
                            rowSums(Zs.deriv * b[id.GK, indRandom, drop = FALSE])
                        else
                            0
                    eta.s <- if (parameterization == "both") 
                        eta.s + c(x$x$WintF.sl[id.GK, , drop = FALSE] %*% Dalpha) * Ys.deriv 
                    else 
                        c(x$x$WintF.sl[id.GK, , drop = FALSE] %*% Dalpha) * Ys.deriv
                }
                ind.K <- rep(unlist(lapply(ind.D, seq_len)), each = nk)
                wk <- unlist(lapply(ind.D, function (n) rep(gaussKronrod(x$control$GKk)$wk, n)))
                P <- c(t(P))
                wkP <- wk * rep(P[!is.na(P)], each = nk)
                Haz[, i] <- exp(eta.tw) * rowsum(xi[ind.K] * wkP * exp(eta.s), id.GK, reorder = FALSE)
            }
            list("survival" = exp(- Haz),
                 "cumulative-Hazard" = Haz,
                 "log-cumulative-Hazard" = log(Haz))
        } else if (method == "spline-PH-GH") {
            T.mat <- matrix(exp(log.survTimes), nrow = n, ncol = nt, byrow = TRUE)
            eta.tw1 <- if (!is.null(W1)) as.vector(W1 %*% gammas) else 0
            gammas.bs <- x$coefficients$gammas.bs
            b <- x$EB$post.b
            wk <- gaussKronrod(x$control$GKk)$wk
            sk <- gaussKronrod(x$control$GKk)$sk
            id.GK <- rep(seq_len(n), each = x$control$GKk)
            Haz <- matrix(0, n, nt)
            for (i in 1:nt) {
                P <- T.mat[, i] / 2
                st <- outer(P, sk + 1)
                data.id <- x$data.id[id.GK, ]
                data.id[[x$timeVar]] <- pmax(c(t(st)) - x$y$lag, 0)
                if (parameterization %in% c("value", "both")) {
                    mfX <- model.frame(x$termsYx, data = data.id)
                    mfZ <- model.frame(x$termsYz, data = data.id)
                    Xs <- model.matrix(x$formYx, mfX)
                    Zs <- model.matrix(x$formYz, mfZ)
                    Ys <- c(Xs %*% x$coefficients$betas) + rowSums(Zs * b[id.GK, , drop = FALSE])
                    eta.s <- c(x$x$WintF.vl[id.GK, , drop = FALSE] %*% alpha) * Ys
                }
                if (parameterization %in% c("slope", "both")) {
                    mfX.deriv <- model.frame(x$termsYx.deriv, data = data.id)
                    mfZ.deriv <- model.frame(x$termsYz.deriv, data = data.id)
                    Xs.deriv <- model.matrix(derivForm$fixed, mfX.deriv)
                    Zs.deriv <- model.matrix(derivForm$random, mfZ.deriv)
                    Ys.deriv <- c(Xs.deriv %*% x$coefficients$betas[indFixed]) +
                        if (length(indRandom) > 1 || indRandom) 
                            rowSums(Zs.deriv * b[id.GK, indRandom, drop = FALSE])
                        else
                            0
                    eta.s <- if (parameterization == "both")
                        eta.s + c(x$x$WintF.sl[id.GK, , drop = FALSE] %*% Dalpha) * Ys.deriv
                    else
                        c(x$x$WintF.sl[id.GK, , drop = FALSE] %*% Dalpha) * Ys.deriv
                }
                strt.s <- rep(x$y$strata, each = x$control$GKk)
                split.Time <- split(c(t(st)), strt.s)
                W2s <- mapply(function (k, t) splineDesign(k, t, ord = x$control$ord, outer.ok = TRUE), 
                    x$control$knots, split.Time, SIMPLIFY = FALSE)
                W2s <- mapply(function (w2s, ind) {
                    out <- matrix(0, n * x$control$GKk, ncol(w2s))
                    out[strt.s == ind, ] <- w2s
                    out
                }, W2s, levels(strt.s), SIMPLIFY = FALSE)
                W2s <- do.call(cbind, W2s)
                eta.ws <- c(W2s %*% gammas.bs)
                Haz[, i] <- exp(eta.tw1) * P * rowsum(wk * exp(eta.ws + eta.s), id.GK, reorder = FALSE)
            }
            list("survival" = exp(- Haz),
                 "cumulative-Hazard" = Haz,
                 "log-cumulative-Hazard" = log(Haz))
        } else {
            W2 <- splineDesign(x$knots, log.survTimes, ord = x$control$ord)
            Y <- c(x$x$Xtime %*% x$coefficients$betas + x$EB$Ztimeb)
            eta <- apply(W2, 1, function (x) {
                w <- matrix(x, n, length(x), TRUE)
                WW <- if (is.null(W1)) w else cbind(w, W1)
                c(WW %*% gammas) + Y * alpha
            })
            list("survival" = exp(- exp(eta)), "cumulative-Hazard" = exp(eta), "log-cumulative-Hazard" = eta)
        }
    }
    if (!return)
        one.fig <- prod(par("mfcol")) == 1
    if (ask && !return) {
        op <- par(ask = TRUE)
        on.exit(par(op))
    }
    if (show[1] && !return) {
        fitY <- fitted(x, process = "Longitudinal", type = "Subject")
        resY <- residuals(x, process = "Longitudinal", type = "Subject")
        plot(fitY, resY, xlab = "Fitted Values", ylab = "Residuals", main = main, ...)
        if (add.smooth) {
            abline(h = 0, lty = 3, col = "grey", lwd = 2)
            panel.smooth(fitY, resY, lwd = 2)
        }
        mtext(caption[1], 3, 0.25, cex = cex.caption)
    } 
    if (show[2] && !return) {
        resY <- residuals(x, process = "Longitudinal", type = "stand-Subject")
        qqnorm(resY, ylab = "Standardized Residuals", main = main, ...)
        if (add.qqline)
            qqline(resY, lty = 3, col = "grey50")
        mtext(caption[2], 3, 0.25, cex = cex.caption)
    }
    if (show[3] && !return) {
        strata <- if (is.null(x$y$strata)) gl(1, n) else x$y$strata
        yy <- rowsum(fitT[["survival"]], strata) / as.vector(table(strata))
        if (add.KM) {
            Time <- exp(x$y$logT)
            failure <- x$y$d
            sf <- survfit(Surv(Time, failure) ~ strata)
            plot(sf, xlab = "Time", ylab = "Survival", main = main, mark.time = FALSE)
            matlines(survTimes, t(yy), ...)
        } else {
            matplot(survTimes, t(yy), xlab = "Time", ylab = "Survival", main = main, ylim = c(0, 1), type = "l", ...)
        }
        mtext(caption[3], 3, 0.25, cex = cex.caption)
    }
    if (show[4] && !return) {
        strata <- if (is.null(x$y$strata)) gl(1, n) else x$y$strata
        yy <- rowsum(fitT[["cumulative-Hazard"]], strata) / as.vector(table(strata))
        matplot(survTimes, t(yy), xlab = "Time", ylab = "Cumulative Hazard", main = main, type = "l", ...)
        mtext(caption[4], 3, 0.25, cex = cex.caption)
    }
    if (show[5] && !return) {
        strata <- if (is.null(x$y$strata)) gl(1, n) else x$y$strata
        yy <- rowsum(fitT[["log-cumulative-Hazard"]], strata) / as.vector(table(strata))
        plot(survTimes, t(yy), xlab = "Time", ylab = "log Cumulative Hazard", main = main, type = "l", ...)
        mtext(caption[5], 3, 0.25, cex = cex.caption)
    }
    if (show[6] && !return) {
        lambda0 <- x$coefficients$lambda0
        plot(lambda0[, "time"], lambda0[, "basehaz"], xlab = "Time", ylab = "", main = main, ...)
        if (add.smooth) {
            panel.smooth(lambda0[, "time"], lambda0[, "basehaz"], lwd = 2)
        }
        mtext(caption[6], 3, 0.25, cex = cex.caption)
    }
    if (show[7] && !return) {
        lambda0 <- x$coefficients$lambda0
        plot(lambda0[, "time"], cumsum(lambda0[, "basehaz"]), xlab = "Time", ylab = "", main = main, type = "s", ...)
        mtext(caption[7], 3, 0.25, cex = cex.caption)
    }
    if (show[8] && !return) {
        yy <- t(fitT[["survival"]])
        matplot(survTimes, yy[, ids], type = "l", col = "black", lty = 1, 
            xlab = "Time", ylab = "Survival", main = main, ylim = c(0, 1), ...)
        mtext(caption[8], 3, 0.25, cex = cex.caption)
    }
    if (show[9] && !return) {
        yy <- t(fitT[["cumulative-Hazard"]])
        matplot(survTimes, yy[, ids], type = "l", col = "black", lty = 1, 
            xlab = "Time", ylab = "Cumulative Hazard", main = main, ...)
        mtext(caption[9], 3, 0.25, cex = cex.caption)
    }
    if (show[10] && !return) {
        yy <- t(fitT[["log-cumulative-Hazard"]])
        matplot(survTimes, yy[, ids], type = "l", col = "black", lty = 1, 
            xlab = "Time", ylab = "log Cumulative Hazard", main = main, ...)
        mtext(caption[10], 3, 0.25, cex = cex.caption)
    }
    if (return)
        invisible(c(fitT, list(survTimes = survTimes)))
    else 
        invisible()
}
