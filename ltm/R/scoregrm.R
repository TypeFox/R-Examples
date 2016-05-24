scoregrm <-
function (thetas, constrained) {
    betas <- betas.grm(thetas, constrained, ind1, ind2, p)
    k <- length(Z)
    cpr <- cprobs(betas, Z)
    diff.cprs <- lapply(cpr, function (x) rbind(x[1, ], diff(x)))
    log.diff.cprs <- lapply(diff.cprs, log)
    log.p.xz <- matrix(0, nfreqs, k)
    for (j in 1:p) {
        log.pr <- log.diff.cprs[[j]]
        xj <- X[, j]
        na.ind <- is.na(xj)
        log.pr <- log.pr[xj, ]
        if (any(na.ind))
            log.pr[na.ind, ] <- 0
        log.p.xz <- log.p.xz + log.pr
    }
    p.xz <- exp(log.p.xz)
    p.x <- c(p.xz %*% GHw)
    p.zx <- p.xz / p.x
    prod.cprs <- lapply(cpr, function (x) x * (1 - x))
    sum.cprs <- lapply(cpr, function (x) {
        nr <- nrow(x)
        1 - rbind(x[1, ], x[-nr, ] + x[-1, ])
    })
    scores.alpha <- lapply(ncatg - 1, numeric)
    scores.beta <- numeric(p)
    scores <- vector("list", p)
    jac <- jacobian(thetas, constrained, ind1, ind2, p)
    for (j in 1:p) {
        pr1 <- diff.cprs[[j]]
        pr2 <- prod.cprs[[j]]
        pr3 <- sum.cprs[[j]]
        for (h in seq(1, ncatg[j] - 1)) {
            mat <- matrix(0, k, nfreqs)
            na.ind <- is.na(X[, j])
            Ind1 <- X[, j] == h
            if (any(na.ind))
                Ind1[na.ind] <- FALSE
            Ind2 <- X[, j] == h + 1
            if (any(na.ind))
                Ind2[na.ind] <- FALSE
            mat[, Ind1] <- pr2[h, ] / pr1[h, ]
            mat[, Ind2] <- -pr2[h, ] / pr1[h + 1, ]
            scores.alpha[[j]][h] <- -sum((p.zx * t(mat) * obs) %*% GHw)
        }
        scores.alpha[[j]] <- scores.alpha[[j]] %*% jac[[j]]
        scores.beta[j] <- sum((p.zx * pr3[X[, j], ] * obs) %*% (Z * GHw), na.rm = TRUE)
        scores[[j]] <- if (constrained) scores.alpha[[j]] else c(scores.alpha[[j]], scores.beta[j])
    }
    if (constrained)
        c(unlist(scores), sum(scores.beta))
    else
        unlist(scores)
}
