flexCPH.fit <-
function (logT, d, X, init.thetas = NULL, control) {
    lgLwph <- function (thetas) {
        gamas <- thetas[1:nk]
        gamas[1:nk] <- cumsum(c(gamas[1], exp(gamas[2:nk])))
        betas <- thetas[-c(1:nk)]
        eta <- as.vector(W %*% gamas + X %*% betas)
        sc <- as.vector(S %*% diff(gamas))
        lgL <- d * (log(sc) + eta - logT) - exp(eta)
        - sum(lgL, na.rm = TRUE)
    }
    Scwph <- function (thetas) {
        gamas <- thetas[1:nk]
        gamas[1:nk] <- cumsum(c(gamas[1], exp(gamas[2:nk])))
        betas <- thetas[-c(1:nk)]
        eta <- as.vector(W %*% gamas + X %*% betas)
        sc <- as.vector(S %*% diff(gamas))
        ew <- - exp(eta)
        out <- (d + ew) * WX
        out[, 1:nk] <- out[, 1:nk] + d * SS / sc
        out <- colSums(out, na.rm = TRUE)
        out[1:nk] <- out[1:nk] %*% jacobian(thetas[1:nk])
        - out
    }
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
    S <- splineDesign(kn[-c(1, length(kn))], logT, ord = control$ord - 1)
    S <- control$ord * S / rep(diff(kn, lag = control$ord + 1), each = length(d))
    ncs <- ncol(S)
    SS <- cbind(- S[, 1], S[, 1:(ncs - 1)] - S[, 2:ncs], S[, ncs])
    nk <- ncol(W)
    WX <- cbind(W, X)
    if (is.null(init.thetas)) {
        init.thetas <-
        init.thetas <- c(-1, seq(-0.5, 0.5, length.out = nk - 1), rep(0, ncol(X)))
    }
    if (is.null(control$parscale))
        control$parscale <- rep(0.01, length(init.thetas))
    opt <- optim(init.thetas, lgLwph, Scwph, method = "BFGS", 
        control = list(maxit = 5000, parscale = control$parscale, reltol = 1e-09))
    thetas <- opt$par
    gamas <- thetas[1:nk]
    gamas[1:nk] <- cumsum(c(gamas[1], exp(gamas[2:nk])))
    betas <- thetas[-(1:nk)]
    eta <- as.vector(W %*% gamas + X %*% betas)
    ew <- exp(eta)
    surv <- exp(- ew)
    CH <- exp(ew)
    logCH <- eta
    H <- if (control$numeriDeriv == "fd") {
        fd.vec(opt$par, Scwph, eps = control$eps.Hes)
    } else {
        cd.vec(opt$par, Scwph, eps = control$eps.Hes)
    }
    if (any(is.na(H) | !is.finite(H))) {
        warning("infinite or missing values in Hessian at convergence.\n")
    } else {
        ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
        if (!all(ev >= -1e-06 * abs(ev[1]))) 
            warning("Hessian matrix at convergence is not positive definite.\n")
    }
    names(gamas) <- paste("bs.", 1:nk, sep = "")
    names(betas) <- colnames(X)
    nams <- c(names(gamas), names(betas))
    dimnames(H) <- list(nams, nams)
    list(coefficients = list(gammas = gamas, betas = betas), Hessian = H, logLik = -opt$value, logT = logT, 
         d = d, X = X, knots = kn, survival = surv, cumHazard = CH, "log.cumHazard" = logCH)
}
