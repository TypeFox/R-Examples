EM <-
function (betas, constraint, iter, verbose = FALSE) {
    for (it in 1:iter) {
        pr <- probs(Z %*% t(betas))
        qr <- 1 - pr
        dvar <- pr * qr
        p.xz <- exp(X %*% t(log(pr)) + mX %*% t(log(qr)))
        p.x <- c(p.xz %*% GHw)
        lgLik <- sum(log(rep(p.x, obs)))
        if (verbose)
            cat("EM iteration:", it, "  -logLik:", -lgLik, "\n")
        p.zx <- (p.xz / p.x) * obs
        Nt <- GHw * colSums(p.zx)
        nb <- matrix(0, p, q.)
        for (i in 1:p) {
            ind <- na.ind[, i]
            Y <- outer(X[, i], pr[, i], "-")
            Y[ind, ] <- 0
            sc <- colSums((p.zx * Y) %*% (Z * GHw))
            hes <- crossprod(Z, (dvar[, i] * Nt) * Z)
            nb[i, ] <- betas[i, ] + solve(hes, sc)
        }
        betas <- nb
        if (!is.null(constraint))
            betas[constraint[, 1:2]] <- constraint[, 3]
    }
    betas
}
