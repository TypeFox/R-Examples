update.logLik.Laplace <-
function (b, betas, sigma, gammas, gammas.bs, alpha, D) {
    log.p.yt <- numeric(n)
    eta.yx <- as.vector(X %*% betas)
    eta.yxT <- as.vector(Xtime %*% betas)
    eta.s <- as.vector(Xs %*% betas)
    eta.tw1 <- as.vector(W1 %*% gammas)
    eta.tw2 <- as.vector(W2 %*% gammas.bs)
    eta.ws <- as.vector(W2s %*% gammas.bs)
    environment(fn) <- environment(gr) <- environment()
    hes.b <- vb <- matrix(0, n, ncz * ncz)
    bb <- b
    for (i in 1:n) {
        # individual values
        ind.i <- id == i
        id.GKi <- id.GK == i
        yi <- y[ind.i]
        eta.yxi <- eta.yx[ind.i]
        eta.yxTi <- eta.yxT[i]
        eta.si <- eta.s[id.GKi]
        Z.ind.i <- Z[ind.i, , drop = FALSE]
        Ztime.i <- Ztime[i, ]
        Zsi <- Zs[id.GKi, , drop = FALSE]
        yi.eta.yxi <- yi - eta.yxi
        logTi <- logT[i]
        eta.tw1i <- eta.tw1[i]
        eta.tw2i <- eta.tw2[i]
        Pi <- P[i]
        eta.wsi <- eta.ws[id.GKi]
        # posterior modes
        opt <- try(optim(b[i, ], fn, gr, method = "BFGS", hessian = TRUE), TRUE)
        if (inherits(opt, "try-error")) {
            opt <- list(par = b[i, ], hessian = matrix(attr(new.b, "hes.b")[i, ], ncz, ncz))
        }
        H <- opt$hessian
        var.b <- if (!inherits(var.b <- try(solve(H), TRUE), "try-error")) var.b else ginv(H)
        bb[i, ] <- opt$par
        vb[i, ] <- c(var.b)
        hes.b[i, ] <- c(H)
        # likelihood contributions
        mu.y.b <- eta.yxi + rowSums(Z.ind.i * rep(bb[i, ], each = ni[i]))
        Yi <- alpha * (eta.si + rowSums(Zsi * rep(bb[i, ], each = GKk)))
        log.p.y.b <- sum(dnorm(yi, mu.y.b, sigma, log = TRUE))
        log.p.t.b <- if (d[i]) {
            eta.tw1i + eta.tw2i + alpha * (eta.yxT[i] + sum(Ztime.i * bb[i, ])) - exp(eta.tw1i) * Pi * sum(wk * exp(eta.wsi + Yi))
        } else {
            - exp(eta.tw1i) * Pi * sum(wk * exp(eta.wsi + Yi))
        }
        log.p.b <- if (diag.D) {
            sum(dnorm(bb[i, ], 0, sqrt(D), log = TRUE))
        } else {
            dmvnorm(bb[i, ], rep(0, ncz), D, log = TRUE)
        }
        log.p.yt[i] <- (log.p.y.b + log.p.t.b + log.p.b) - 0.5 * log(det(H))
    }
    res <- cons.logLik + sum(log.p.yt, na.rm = TRUE)
    attr(res, "b") <- bb
    attr(res, "vb") <- vb
    attr(res, "hes.b") <- hes.b
    res
}
