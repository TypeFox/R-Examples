vcov.weibull.frailty <-
function (object, sand.se = FALSE, ...) {
    inv.H <- solve(object$hessian)
    betas <- object$coefficients$betas
    shape <- object$coefficients$shape
    scale <- object$coefficients$scale
    var.fr <- object$coefficients$var.frailty
    out <- if (sand.se) {
        logT <- log(object$y[, 1])
        d <- object$y[, 2]
        id <- object$id
        Time <- exp(logT)
        X <- object$x
        p <- ncol(X)
        N <- nrow(X)
        n <- length(unique(id))
        D <- as.vector(tapply(d, id, sum))
        if (!ncol(X))
            X <- as.matrix(rep(0, N))
        Xd <- rowsum(X * d, id, FALSE)
        sum.d <- as.vector(tapply(d, id, sum))
        sum.dlogT <- as.vector(tapply(d * logT, id, sum))
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
        sc.betas <- (Xd - (D + theta) * X.Lambda0.eta / theta.Lambda0.eta)
        sc.scale <- scale * c(sum.d / scale - (D + theta) * P[, 2] / theta.Lambda0.eta)
        sc.shape <- shape * (sum.d / shape + sum.dlogT - (D + theta) * P[, 3] / theta.Lambda0.eta)
        sc.var.fr <- - theta * (log(theta) + 1 - digamma(theta) + digamma(D + theta) - 
            log.theta.Lambda0.eta - (D + theta) / theta.Lambda0.eta)
        score <- if (p > 0) cbind(sc.betas, sc.scale, sc.shape, sc.var.fr) else cbind(sc.scale, sc.shape, sc.var.fr)
        out.score <- colSums(t(apply(score, 1, function (x) x %o% x)))
        dim(out.score) <- dim(inv.H)
        out.score <- 0.5 * (out.score + t(out.score))
        inv.H %*% out.score %*% inv.H
    } else
        inv.H
    out <- 0.5 * (out + t(out))
    nams <- c(names(betas), "Log(scale)", "Log(shape)", "Log(var.frlty)")
    dimnames(out) <- list(nams, nams)
    out
}
