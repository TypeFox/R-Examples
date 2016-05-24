initial.survOld <-
function (logT, d, X, parameterization, method, control, extra = NULL) {
    old <- options(warn = (-1))
    on.exit(options(old))
    if (method == "Cox-PH-GH") {
        Time <- exp(logT)
        cph <- if (parameterization == "value") {
            start <- extra$times
            DD <- data.frame(id = extra$id, start = start, Time = Time[extra$id])
            stop <- unlist(lapply(split(DD, DD$id), function (x) c(x$start[-1], x$Time[1])))
            event <- ave(d[extra$id], extra$id, FUN = function(x) {
                nx <- length(x)
                if (nx == 1) x[1] else c(rep(0, nx-1), x[1]) 
            })
            WW <- if ((nx <- ncol(X)) == 1) extra$y else cbind(X[extra$id, -nx], extra$y)
            coxph(Surv(start, stop, event) ~ WW)
        } else {
            coxph(Surv(Time, d) ~ X)
        }
        coefs <- coef(cph)
        nk <- length(coefs)
        out <- switch(parameterization, 
            "value" = if (nk > 1) list(gammas = coefs[-nk], alpha = coefs[nk]) else list(alpha = coefs[nk]),
            "slope" = if (nk > 1) list(gammas = coefs[-nk], Dalpha = coefs[nk]) else list(Dalpha = coefs[nk]),
            "both" = if (nk > 2) 
                list(gammas = coefs[-c(nk-1, nk)], alpha = coefs[nk-1], Dalpha = coefs[nk]) else list(alpha = coefs[nk-1], Dalpha = coefs[nk])
        )
        lambda0 <- basehaz(cph, FALSE)
        out <- c(out, list(lambda0 = lambda0$hazard))
        out
    } else if (method == "weibull-PH-GH") {
        dat <- data.frame(Time = exp(logT), d = d)
        init.fit <- survreg(Surv(Time, d) ~ X, data = dat)
        coefs <- - init.fit$coef / init.fit$scale
        nk <- length(coefs)
        switch(parameterization, 
            "value" = list(gammas = coefs[-nk], alpha = coefs[nk], sigma.t = 1 / init.fit$scale),
            "slope" = list(gammas = coefs[-nk], Dalpha = coefs[nk], sigma.t = 1 / init.fit$scale),
            "both" = list(gammas = coefs[1:(nk-2)], alpha = coefs[nk-1], Dalpha = coefs[nk], sigma.t = 1 / init.fit$scale)
        )
    } else if (method == "weibull-AFT-GH") {
        dat <- data.frame(Time = exp(logT), d = d)
        init.fit <- survreg(Surv(Time, d) ~ X, data = dat)
        coefs <- - init.fit$coef
        nk <- length(coefs)
        switch(parameterization, 
            "value" = list(gammas = coefs[-nk], alpha = coefs[nk], sigma.t = 1 / init.fit$scale),
            "slope" = list(gammas = coefs[-nk], Dalpha = coefs[nk], sigma.t = 1 / init.fit$scale),
            "both" = list(gammas = coefs[1:(nk-2)], alpha = coefs[nk-1], Dalpha = coefs[nk], sigma.t = 1 / init.fit$scale)
        )
    } else if (method == "piecewise-PH-GH") {
        n <- nrow(d)
        p <- ncol(d)
        XX <- X[rep(seq_len(n), p), , drop = FALSE]; rownames(XX) <- NULL
        dat <- data.frame(Time = c(logT), d = c(d), xi = gl(p, n), XX)
        init.fit <- glm(d ~ . + offset(log(Time)) - Time - 1, family = poisson, data = dat[dat$Time > 0, ])
        out <- list(xi = exp(init.fit$coefficients[seq_len(length(control$knots) + 1)]))
        coefs <- if (parameterization == "value") {
            start <- extra$times
            DD <- data.frame(id = extra$id, start = start, Time = extra$Time[extra$id])
            stop <- unlist(lapply(split(DD, DD$id), function (x) c(x$start[-1], x$Time[1])))
            dd <- extra$d
            event <- ave(dd[extra$id], extra$id, FUN = function (x) {
                nx <- length(x)
                if (nx == 1) x[1] else c(rep(0, nx-1), x[1]) 
            })
            WW <- if ((nx <- ncol(X)) == 1) extra$y else cbind(X[extra$id, -nx], extra$y)
            coxph(Surv(start, stop, event) ~ WW)$coefficients
        } else {
            coxph(Surv(Time, d) ~ X, as.data.frame(extra[c("Time", "d")]))$coefficients
        }
        nk <- length(coefs)
        out <- c(out, switch(parameterization, 
            "value" = if (nk > 1) list(gammas = coefs[-nk], alpha = coefs[nk]) else list(alpha = coefs[nk]),
            "slope" = if (nk > 1) list(gammas = coefs[-nk], Dalpha = coefs[nk]) else list(Dalpha = coefs[nk]),
            "both" = if (nk > 2) 
                list(gammas = coefs[1:(nk-2)], alpha = coefs[nk-1], Dalpha = coefs[nk]) else list(alpha = coefs[nk-1], Dalpha = coefs[nk])
        ))
        out
    } else if (method == "spline-PH-GH" || method == "spline-PH-Laplace") {
        dat <- data.frame(Time = exp(logT), d = d)
        init.fit <- survreg(Surv(Time, d) ~ X, data = dat)
        coefs <- - init.fit$coef / init.fit$scale
        nk <- length(coefs)
        out <- switch(parameterization, 
            "value" = if (nk > 2) list(gammas = coefs[-c(1, nk)], alpha = coefs[nk]) else list(alpha = coefs[nk]),
            "slope" = if (nk > 2) list(gammas = coefs[-c(1, nk)], Dalpha = coefs[nk]) else list(Dalpha = coefs[nk]),
            "both" = if (nk > 3) 
                list(gammas = coefs[-c(1, nk-1, nk)], alpha = coefs[nk-1], Dalpha = coefs[nk]) else list(alpha = coefs[nk-1], Dalpha = coefs[nk])
        )
        xi <- 1 / init.fit$scale
        phi <- exp(coefs[1])
        logh <- log(phi * xi * exp(logT)^(xi - 1))
        out$gammas.bs <- as.vector(lm.fit(extra$W2, logh)$coefficients)
        out
    } else if (method == "ch-Laplace" || method == "ch-GH") {
        dat <- data.frame(Time = exp(logT), d = d)
        init.fit <- survreg(Surv(Time, d) ~ X, data = dat)
        coefs <- - coef(init.fit) / init.fit$scale
        min.x <- min(logT)
        max.x <- max(logT)
        kn <- if (is.null(control$knots)) {
            kk <- seq(0, 1, length.out = control$lng.in.kn + 2)[-c(1, control$lng.in.kn + 2)]
            quantile(logT[d == 1], kk, names = FALSE)
        } else {
            control$knots
        }
        kn <- sort(c(rep(c(min.x, max.x), control$ord), kn))
        W <- splineDesign(kn, logT, ord = control$ord)
        nk <- ncol(W)
        nx <- NCOL(X)
        logH <- coefs[1] + logT / init.fit$scale
        coefs <- c(as.vector(lm.fit(W, logH)$coefficients), coefs[-1])
        list(gammas = c(sort(coefs[1:nk]), if (nx > 1) coefs[seq(nk + 1, nk + nx - 1)] else NULL), 
            alpha = coefs[nk + nx])
    }
}
