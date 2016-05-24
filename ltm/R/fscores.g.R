fscores.g <-
function (betas, X, method) {
    logf.z <- function (z, y, betas) {
        gammas <- lapply(betas, function (x) {
            n <- length(z)
            nx <- length(x)
            c(plogis(matrix(x[-nx], n, nx - 1, TRUE) - x[nx] * z), 1)
        })
        log.prs <- lapply(gammas, function (x) {
            nc <- length(x)
            log(c(x[1], x[2:nc] - x[1:(nc - 1)]))
        })
        log.pxz <- numeric(p)
        for (i in 1:p) {
            log.pxz[i] <- if (!is.na(y[i])) log.prs[[i]][y[i]] else 0
        }
        if (prior)
            - (sum(log.pxz, na.rm = TRUE) + dnorm(z, log = TRUE)) 
        else
            - sum(log.pxz, na.rm = TRUE)
    }
    fscore <- function (logf.z, y, betas) {
        opt <- optim(0.0, fn = logf.z, method = "BFGS", hessian = TRUE, y = y, betas = betas)
        hc <- c(1/opt$hes)
        list(mu = opt$par, hes = hc)
    }
    if (method == "EB") {
        scores.ML <- hes.ML <- numeric(nx)
        for (i in 1:nx) {
            out <- fscore(logf.z = logf.z, y = X[i, ], betas = betas)
            scores.ML[i] <- out$mu
            hes.ML[i] <- out$hes
        }
        res$z1 <- scores.ML
        res$se.z1 <- sqrt(hes.ML)
    }
    if (method == "EAP") {
        Z <- object$GH$Z
        GHw <- object$GH$GHw
        cpr <- cprobs(betas, Z)
        diff.cprs <- lapply(cpr, function (x) rbind(x[1, ], diff(x)))
        log.diff.cprs <- lapply(diff.cprs, log)
        log.p.xz <- matrix(0, nrow(X), length(Z))
        for (j in 1:p) {
            log.pr <- log.diff.cprs[[j]]
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
        res$z1 <- c(p.zx %*% (Z * GHw))
        res$se.z1 <- sqrt(c(p.zx %*% (Z * Z * GHw)) - res$z1^2)        
    }
    if (method == "MI") {
        constrained <- object$constrained
        p <- length(betas)
        ncatg <- sapply(betas, length)
        if (constrained)
            betas[-p] <- lapply(betas[-p], function (x) x[-length(x)])
        betas <- unlist(betas, use.names = FALSE)
        ind1 <- if (constrained) c(1, cumsum(ncatg[-p] - 1) + 1) else c(1, cumsum(ncatg[-p]) + 1)
        ind2 <- if (constrained) cumsum(ncatg - 1) else cumsum(ncatg)        
        Var.betas <- solve(object$hessian)
        scores.B <- hes.B <- array(0, dim = c(nx, B))
        for (b in 1:B) {
            betas. <- mvrnorm(1, betas, Var.betas)
            betas. <- betas.grm(betas., constrained, ind1, ind2, p, trasform = FALSE)
            for (i in 1:nx) {
                out <- fscore(logf.z = logf.z, y = X[i, ], betas = betas.)
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
