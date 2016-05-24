fscores.r <-
function (betas, X, method) {
    logf.z <- function (z, y, betas) {
        pr <- probs(c(betas %*% c(1, z)))
        if (prior)
            -sum(dbinom(y, 1, pr, log = TRUE), na.rm = TRUE) - dnorm(z, log = TRUE)
        else
            -sum(dbinom(y, 1, pr, log = TRUE), na.rm = TRUE)
    }
    fscore <- function (logf.z, y, betas) {
        opt <- optim(0, fn = logf.z, method = "BFGS", hessian = TRUE, y = y, betas = betas)
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
        pr <- probs(Z %*% t(betas))
        mX <- 1 - X
        if (any(na.ind <- is.na(X)))
            X[na.ind] <- mX[na.ind] <- 0
        p.xz <- exp(X %*% t(log(pr)) + mX %*% t(log(1 - pr)))
        p.x <- c(p.xz %*% GHw)
        p.zx <- p.xz / p.x
        res$z1 <- c(p.zx %*% (Z[, 2] * GHw))
        res$se.z1 <- sqrt(c(p.zx %*% (Z[, 2] * Z[, 2] * GHw)) - res$z1^2)
    }
    if (method == "MI") {
        constraint <- object$constraint
        betas <- c(betas[, 1], betas[1, 2])
        if (!is.null(constraint))
            betas <- betas[-constraint[, 1]]
        Var.betas <- vcov(object, robust.se)
        scores.B <- hes.B <- array(0, dim = c(nx, B))
        for (b in 1:B) {
            betas. <- mvrnorm(1, betas, Var.betas)
            betas. <- betas.rasch(betas., constraint, p)
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
        SV <- rowSums(SV)/(B - 1)
        hes.av <- hes.av + (1 + 1/B) * SV
        res$z1 <- scores.av
        res$se.z1 <- sqrt(hes.av)
        attr(res, "zvalues.MI") <- scores.B
        attr(res, "var.zvalues.MI") <- hes.B
    }
    res
}
