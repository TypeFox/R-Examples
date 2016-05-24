fscores.gp <-
function (betas, X, method) {
    logf.z <- function (z, y, betas) {
        log.prs <- crf.GPCM(betas, z, IRT.param = object$IRT.param, log = TRUE)
        log.pxz <- numeric(p)
        for (i in 1:p) {
            log.pxz[i] <- if (!is.na(y[i])) log.prs[[i]][y[i]] else 0
        }
        if (prior)
            - (sum(log.pxz, na.rm = TRUE) + dnorm(z, log = TRUE))
        else
            - sum(log.pxz, na.rm = TRUE)
    }
    fscore <- function (logf.z, y, betas, start) {
        opt <- optim(start, fn = logf.z, method = "BFGS", hessian = TRUE, y = y, betas = betas)
        hc <- c(1/opt$hes)
        list(mu = opt$par, hes = hc)
    }
    Z <- object$GH$Z
    GHw <- object$GH$GHw
    log.crf <- crf.GPCM(betas, Z, object$IRT.param, log = TRUE)
    log.p.xz <- matrix(0, nrow(X), length(Z))
    for (j in 1:p) {
        log.pr <- log.crf[[j]]
        xj <- X[, j]
        na.ind <- is.na(xj)
        log.pr <- log.pr[xj, , drop = FALSE]
        if (any(na.ind))
            log.pr[na.ind, ] <- 0
        log.p.xz <- log.p.xz + log.pr
    }
    p.xz <- exp(log.p.xz)
    p.x <- c(p.xz %*% GHw)
    p.zx <- p.xz / p.x
    if (method == "EB") {
        z.st <- c(p.zx %*% (Z * GHw))
        scores.ML <- hes.ML <- numeric(nx)
        for (i in 1:nx) {
            out <- fscore(logf.z = logf.z, y = X[i, ], betas = betas, start = z.st[i])
            scores.ML[i] <- out$mu
            hes.ML[i] <- out$hes
        }
        res$z1 <- scores.ML
        res$se.z1 <- sqrt(hes.ML)
    }
    if (method == "EAP") {
        res$z1 <- c(p.zx %*% (Z * GHw))
        res$se.z1 <- sqrt(c(p.zx %*% (Z * Z * GHw)) - res$z1^2)        
    }
    if (method == "MI") {
        constraint <- object$constraint
        p <- length(betas)
        ncatg <- sapply(betas, length)
        vec.betas <- if (constraint == "gpcm") {
            unlist(betas, use.names = FALSE)
        } else if (constraint == "1PL") {
            betas[seq(1, p - 1)] <- lapply(betas[seq(1, p - 1)], function (x) x[-length(x)])
            unlist(betas, use.names = FALSE)
        } else {
            betas <- lapply(betas, function (x) x[-length(x)])
            unlist(betas, use.names = FALSE)
        }
        Var.betas <- vcov(object, robust.se =  robust.se)
        z.st <- c(p.zx %*% (Z * GHw))
        scores.B <- hes.B <- array(0, dim = c(nx, B))
        for (b in 1:B) {
            new.betas <- mvrnorm(1, vec.betas, Var.betas)
            new.betas <- betas.gpcm(new.betas, p, ncatg, constraint)
            for (i in 1:nx) {
                out <- fscore(logf.z = logf.z, y = X[i, ], betas = new.betas, start = z.st[i])
                scores.B[i, b] <- out$mu
                hes.B[i, b] <- out$hes
            }
        }
        scores.av <- rowMeans(scores.B)
        hes.av <- rowMeans(hes.B)
        SV <- array(0, dim = c(nx, B))
        for (b in 1:B) {
            for (i in 1:nx) {
                sc.dif <- scores.B[i, b] - scores.av[i]
                SV[i, b] <- outer(sc.dif, sc.dif)
            }
        }
        SV <- rowSums(SV) / (B - 1)
        hes.av <- hes.av + (1 + 1/B) * SV
        res$z1 <- scores.av
        res$se.z1 <- sqrt(hes.av)
        attr(res, "zvalues.MI") <- scores.B
        attr(res, "var.zvalues.MI") <- hes.B
    }
    res
}
