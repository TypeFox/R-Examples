scoregpcmSNW <-
function (object) {
    X <- if (!is.data.frame(object$X)) as.data.frame(object$X) else object$X
    X[] <- lapply(X, factor)
    ncatg <- as.vector(sapply(X, function (x) length(levels(x))))
    X <- data.matrix(X)
    Z <- object$GH$Z
    GHw <- object$GH$GHw
    IRT.param <- object$IRT.param
    constraint <- object$constraint
    n <- nrow(X)
    p <- ncol(X)
    XX <- lapply(1:p, function (j) outer(X[, j], seq(1, ncatg[j] - 1), ">") * 1)
    betas <- object$coefficients
    eta <- eta1 <- log.crf <- num <- den <- vector("list", p)
    log.p.xz <- matrix(0, nrow(X), length(Z))
    na.ind <- is.na(X)
    for (j in 1:p) {
        bt <- betas[[j]]
        nbt <- length(bt)
        eta1[[j]] <- if (IRT.param) t(outer(Z, bt[-nbt], "-")) else matrix(Z, nbt - 1, length(Z), TRUE)
        eta[[j]] <- if (IRT.param) bt[nbt] * eta1[[j]] else outer(bt[-nbt], bt[nbt] * Z , "+")
        num[[j]] <- exp(apply(eta[[j]], 2, cumsum))
        if (!is.matrix(num[[j]]))
            num[[j]] <- t(num[[j]])
        den[[j]] <- 1 + colSums(num[[j]])
        crf <- rbind(1/den[[j]], num[[j]] / rep(den[[j]], each = ncatg[j] - 1))
        eps <- .Machine$double.eps^(1/2)
        if (any(ind <- crf == 1))
            crf[ind] <- 1 - eps
        if (any(ind <- crf == 0))
            crf[ind] <- eps
        log.pr <- log(crf)
        xj <- X[, j]
        log.pr <- log.pr[xj, ]
        if (any(na.ind[, j]))
            log.pr[na.ind[, j], ] <- 0
        log.p.xz <- log.p.xz + log.pr
    }
    p.xz <- exp(log.p.xz)
    p.x <- c(p.xz %*% GHw)
    p.zx <- p.xz / p.x
    check.alpha <- constraint == "gpcm" || constraint == "1PL"
    if (check.alpha)
        scores.alpha <- matrix(0, n, p)
    scores.deltas <- lapply(ncatg - 1, numeric)
    for (j in 1:p) {
        if (check.alpha) {
            xj <- X[, j]
            etaj <- eta[[j]]
            eta1j <- apply(eta1[[j]], 2, cumsum)
            I1 <- rbind(0, eta1j) - rep(colSums(num[[j]] * eta1j) / den[[j]], each = ncatg[j])
            scores.alpha[, j] <- - (p.zx * I1[xj, ]) %*% GHw
        }
        alpha <- betas[[j]][ncatg[j]]
        ii <- seq(1, ncatg[j] - 1)
        ind1 <- unlist(sapply(ii, ":", b = ncatg[j] - 1)) 
        ind2 <- rep(ii, rev(ii))
        Part1 <- if (IRT.param) - alpha * XX[[j]] else XX[[j]]
        Part2 <- if (IRT.param) - alpha * rowsum(num[[j]][ind1, , drop = FALSE], ind2) else rowsum(num[[j]][ind1, , drop = FALSE], ind2)
        scores.deltas[[j]] <- do.call("cbind", lapply(ii, function (i) {
            I2 <- outer(Part1[, i], Part2[i, ] / den[[j]], "-")
            if (any(na.ind[, j]))
                I2[na.ind[, j], ] <- 0
            - (p.zx * I2) %*% GHw
        }))
    }
    score <- if (constraint == "gpcm") {
        do.call("cbind",
            mapply(function (x, y) cbind(x, y), scores.deltas, as.data.frame(scores.alpha), SIMPLIFY = FALSE, USE.NAMES = FALSE))
    } else if (constraint == "1PL") {
        cbind(do.call("cbind", scores.deltas), rowSums(scores.alpha))
    } else {
        do.call("cbind", scores.deltas)
    }
    out.score <- colSums(t(apply(score, 1, function(x) x %o% x)), na.rm = TRUE)
    dim(out.score) <- c(ncol(score), ncol(score))
    0.5 * (out.score + t(out.score))
}
