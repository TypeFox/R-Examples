weibull.frailty.fit <-
function (logT, d, X, id, init.thetas, control = list()) {
    Time <- exp(logT)
    p <- ncol(X)
    N <- nrow(X)
    n <- length(unique(id))
    D <- as.vector(tapply(d, id, sum))
    if (!ncol(X))
        X <- as.matrix(rep(0, N))
    Xd <- colSums(X * d)
    sum.d <- sum(d)
    sum.dlogT <- sum(d * logT)
    con <- list(optimizer = "optim", parscale = NULL, maxit = 500, numeriDeriv = "cd", eps.Hes = 1e-03)
    con[names(control)] <- control
    clnams <- colnames(X)
    dimnames(X) <- names(logT) <- names(Time) <- names(d) <- names(id) <- names(Xd) <- NULL 
    fn <- function (thetas) {
        betas <- thetas[seq_len(p)]
        scale <- exp(thetas[p + 1])
        shape <- exp(thetas[p + 2])
        var.fr <- exp(thetas[p + 3])
        theta <- 1 / var.fr
        eta <- if (p > 0) c(X %*% betas) else rep(0, N)
        log.lambda0 <- log(shape * scale) + (shape - 1) * logT
        Lambda0 <- scale * Time^shape
        P1 <- sum(d * (log.lambda0 + eta))
        P2 <- n * (theta * log(theta) - lgamma(theta))
        P3 <- sum(lgamma(D + theta) - (D + theta) * log(theta + c(tapply(Lambda0 * exp(eta), id, sum))))
        - (P1 + P2 + P3)
    }
    gr <- function (thetas) {
        betas <- thetas[seq_len(p)]
        scale <- exp(thetas[p + 1])
        shape <- exp(thetas[p + 2])
        var.fr <- exp(thetas[p + 3])
        theta <- 1 / var.fr
        eta <- if (p > 0) c(X %*% betas) else rep(0, N)
        exp.eta <- exp(eta)
        log.lambda0 <- log(shape * scale) + (shape - 1) * logT
        Lambda0 <- scale * Time^shape
        Lambda0.eta <- Lambda0 * exp(eta)
        mat.id <- cbind(Lambda0.eta, Time^shape * exp(eta), logT * Lambda0.eta, Lambda0.eta * X)
        P <- rowsum(mat.id, id, FALSE)
        Lambda0.eta <- P[, 1]
        theta.Lambda0.eta <- theta + Lambda0.eta
        log.theta.Lambda0.eta <- log(theta.Lambda0.eta)
        X.Lambda0.eta <- P[, seq(4, ncol(P)), drop = FALSE]
        sc.betas <- - c(Xd - colSums((D + theta) * X.Lambda0.eta / theta.Lambda0.eta))
        sc.scale <- - scale * c(sum.d / scale - sum((D + theta) * P[, 2] / theta.Lambda0.eta))
        sc.shape <- - shape * (sum.d / shape + sum.dlogT - sum((D + theta) * P[, 3] / theta.Lambda0.eta))
        sc.var.fr <- theta * sum(log(theta) + 1 - digamma(theta) + digamma(D + theta) - 
            log.theta.Lambda0.eta - (D + theta) / theta.Lambda0.eta)
        if (p > 0) c(sc.betas, sc.scale, sc.shape, sc.var.fr) else c(sc.scale, sc.shape, sc.var.fr)
    }
    if (is.null(init.thetas) || length(init.thetas) != p + 3)
        init.thetas <- rep(0.01, p + 3)
    names(init.thetas) <- NULL
    psc <- if (is.null(xx <- con$parscale)) rep(c(1, 0.1), c(p, 3)) else xx
    opt <- if (con$optimizer == "optim") {
        optim(init.thetas, fn, gr, method = "BFGS", control = list(maxit = con$maxit, parscale = psc))
    } else {
        nlminb(init.thetas, fn, gr, scale = psc, control = list(iter.max = con$maxit))
    }
    H <- if (con$numeriDeriv == "cd") cd.vec(opt$par, gr, eps = con$eps.Hes) else fd.vec(opt$par, gr, eps = con$eps.Hes)
    if (any(is.na(H) | !is.finite(H))) {
        warning("infinite or missing values in Hessian at convergence.\n")
    } else {
        ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
        if (!all(ev >= -1e-06 * abs(ev[1]))) 
            warning("Hessian matrix at convergence is not positive definite.\n")
    }
    betas <- opt$par[seq_len(p)]
    names(betas) <- clnams
    scale <- exp(opt$par[p + 1])
    shape <- exp(opt$par[p + 2])
    var.fr <- exp(opt$par[p + 3])
    list(coefficients = list(betas = betas, scale = scale, shape = shape, var.frailty = var.fr), hessian = H, 
         logLik = -opt[[2]], control = con)
}
