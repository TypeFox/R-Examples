predict.b <-
function (method, y, X, Xtime, Z, Ztime, betas, sigma, Time, W1, gammas, alpha, sigma.t, D, 
    id, control, knots) {
    WW <- if (method == "Cox-PH-GH") {
        NA
    } else if (method == "weibull-GH") {
        if (is.null(W1)) as.matrix(rep(1, length(Time))) else cbind(1, W1)
    } else {
        logT <- log(Time)
        W2 <- splineDesign(knots, logT, ord = control$ord)
        S <- splineDesign(knots[-c(1, length(knots))], logT, ord = control$ord - 1)
        S <- control$ord * S / rep(diff(knots, lag = control$ord + 1), each = length(Time))
        ncs <- ncol(S)
        #SS <- cbind(- S[, 1], S[, 1:(ncs - 1)] - S[, 2:ncs], S[, ncs])
        nk <- ncol(W2)
        if (is.null(W1)) cbind(W2) else cbind(W2, W1)
    }
    ncx <- ncol(X)
    ncz <- ncol(Z)
    ncww <- ncol(WW)
    n <- length(Time)
    N <- length(y)
    ni <- as.vector(tapply(id, id, length))
    diag.D <- ncz != ncol(D)
    out <- if (method %in% c("Cox-PH-GH", "weibull-GH", "ch-GH", "spline-GH-PH")) {
        GH <- gauher(control$GHk)
        b <- as.matrix(expand.grid(lapply(1:ncz, function (k, u) u$x, u = GH)))
        k <- nrow(b)
        wGH <- as.matrix(expand.grid(lapply(1:ncz, function (k, u) u$w, u = GH)))
        wGH <- 2^(ncz/2) * apply(wGH, 1, prod) * exp(rowSums(b * b)) * control$det.inv.chol.VC
        b <- sqrt(2) * t(control$inv.chol.VC %*% t(b))
        b2 <- if (ncz == 1) b * b else t(apply(b, 1, function (x) x %o% x))
        Ztb <- Z %*% t(b)
        Ztime.b <- Ztime %*% t(b)
        eta.yx <- as.vector(X %*% betas)
        eta.yxT <- as.vector(Xtime %*% betas)
        eta.tw <- as.vector(WW %*% gammas)
        Y <- eta.yxT + Ztime.b
        eta.t <- eta.tw + alpha * Y
        mu.y <- eta.yx + Ztb
        logNorm <- dnorm(y, mu.y, sigma, TRUE)
        log.p.yb <- rowsum(logNorm, id)
        log.p.tb <- if (method == "ph-GH") {
            NA
        } else if (method == "weibull-GH") {
            w <- (log(Time) - eta.t) / sigma.t
            - exp(w)
        } else {
            - exp(eta.t)
        }
        log.p.b <- if (ncz == 1) {
            dnorm(b, sd = sqrt(D), log = TRUE)
        } else {
            if (diag.D) {
                rowSums(dnorm(b, sd = rep(sqrt(D), each = k), log = TRUE))
            } else {
                dmvnorm(b, rep(0, ncz), D, TRUE)
            }
        }
        p.ytb <- exp((log.p.yb + log.p.tb) + rep(log.p.b, each = n))
        dimnames(p.ytb) <- NULL
        p.yt <- c(p.ytb %*% wGH)
        p.byt <- p.ytb / p.yt
        p.byt %*% (b * wGH)
    } else {
        environment(update.bCH) <- environment(fn.b) <- environment(gr.b) <- environment()
        vb <- matrix(0, n, ncz * ncz); cons.logLik <- 0
        new.b <- update.bCH(matrix(0, n, ncz), vb, betas, sigma, c(gammas, alpha), D)
        attr(new.b, "b")
    }
    attr(out, "WW") <- WW
    out
}
