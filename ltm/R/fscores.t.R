fscores.t <-
function (thetas, X, method) {
    logf.z <- function (z, y, thetas) {
        betas <- thetas[, 2:3]
        cs <- plogis(thetas[, 1]) * object$max.guessing
        pr <- cs + (1 - cs) * probs(c(betas %*% c(1, z)))
        if (prior)
            -sum(dbinom(y, 1, pr, log = TRUE), na.rm = TRUE) - dnorm(z, log = TRUE)
        else 
            -sum(dbinom(y, 1, pr, log = TRUE), na.rm = TRUE)
    }
    fscore <- function (logf.z, y, thetas) {
        opt <- optim(0, fn = logf.z, method = "BFGS", hessian = TRUE, y = y, thetas = thetas)
        hc <- c(1 / opt$hes)
        list(mu = opt$par, hes = hc)
    }
    if (method == "EB") {
        scores.ML <- hes.ML <- numeric(nx)
        for (i in 1:nx) {
            out <- fscore(logf.z = logf.z, y = X[i, ], thetas = thetas)
            scores.ML[i] <- out$mu
            hes.ML[i] <- out$hes
        }
        res$z1 <- scores.ML
        res$se.z1 <- sqrt(hes.ML)
    }
    if (method == "EAP") {
        Z <- object$GH$Z
        GHw <- object$GH$GHw
        betas <- thetas[, 2:3]
        cs <- plogis(thetas[, 1]) * object$max.guessing
        cs.mat <- matrix(cs, nrow(Z), nrow(betas), TRUE)
        pr <- cs + (1 - cs) * probs(Z %*% t(betas))
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
        if (!is.null(constraint))
            thetas <- thetas[-((constraint[, 2] - 1) * p + constraint[, 1])]
        thetas <- unique(c(thetas))
        Var.betas <- solve(object$hessian)
        scores.B <- hes.B <- array(0, dim = c(nx, B))
        for (b in 1:B) {
            thetas. <- mvrnorm(1, thetas, Var.betas)
            thetas. <- thetas.tpm(thetas., object$type, constraint, p)
            for (i in 1:nx) {
                out <- fscore(logf.z = logf.z, y = X[i, ], thetas = thetas.)
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
