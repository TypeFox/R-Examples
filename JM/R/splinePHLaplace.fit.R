splinePHLaplace.fit <-
function (x, y, id, initial.values, initial.EB = NULL, control) {
    # response vectors
    logT <- as.vector(y$logT)
    d <- as.vector(y$d)
    y <- as.vector(y$y)
    # design matrices
    X <- x$X
    Xtime <- x$Xtime
    Z <- x$Z
    Ztime <- x$Ztime
    W1 <- x$W
    W2 <- x$W2
    W2s <- x$W2s
    dimnames(X) <- dimnames(Xtime) <- dimnames(Z) <- dimnames(Ztime) <- dimnames(W1) <- NULL
    attr(X, "assign") <- attr(X, "contrasts") <- attr(Xtime, "assign") <- attr(Xtime, "contrasts") <- NULL
    attr(Z, "assign") <- attr(Ztime, "assign") <- NULL
    # sample size settings
    ncx <- ncol(X)
    ncz <- ncol(Z)
    ncw1 <- if (is.null(W2)) 0 else ncol(W1)
    nk <- ncol(W2)
    n <- length(logT)
    N <- length(y)
    ni <- as.vector(tapply(id, id, length))
    # Gauss-Kronrod rule
    P <- as.vector(x$P)
    id.GK <- rep(seq_along(logT), each = control$GKk)
    GKk <- control$GKk
    # crossproducts and others
    XtX <- crossprod(X)
    ZtZ <- lapply(split(Z, id), function (x) crossprod(matrix(x, ncol = ncz))); names(ZtZ) <- NULL
    ZtZ <- matrix(unlist(ZtZ), n, ncz * ncz, TRUE)
    outer.Ztime <- lapply(1:n, function (x) Ztime[x, ] %o% Ztime[x, ])
    # initial values
    betas <- as.vector(initial.values$betas)
    sigma <- initial.values$sigma
    gammas <- as.vector(initial.values$gammas)
    gammas.bs <- as.vector(initial.values$gammas.bs)
    alpha <- as.vector(initial.values$alpha)
    D <- initial.values$D
    diag.D <- !is.matrix(D)
    if (!diag.D) dimnames(D) <- NULL else names(D) <- NULL
    # initialize Laplace approximation components
    b <- trc.y1 <- if (is.null(initial.EB) || !all(dim(initial.EB) == c(n, ncz))) matrix(0, n, ncz) else initial.EB
    dimnames(b) <- NULL
    trc.y2 <- trc.y3 <- matrix(0, n, ncz * ncz, TRUE)
    trc.t1 <- trc.t2 <- numeric(n)
    cons.logLik <- 0.5 * n * ncz * log(2 * pi)
    # Fix environments for functions
    #environment(update.logLik.Laplace) <- environment(fn) <- environment(gr) <- environment()
    #environment(logsurvCH) <- environment(ScsurvCH) <- environment(SclongCH) <- environment()
    #environment(LogLik.chLaplace) <- environment(Score.chLaplace) <- environment()
    #old <- options(warn = (-1))
    #on.exit(options(old))    
    # EM iterations
    iter <- control$iter.EM
    Y.mat <- matrix(0, iter, ncx + 1)
    T.mat <- matrix(0, iter, ncw1 + nk + 1)
    B.mat <- if (diag.D) matrix(0, iter, ncz) else matrix(0, iter, ncz * ncz)
    lgLik <- numeric(iter)
    conv <- FALSE
    new.b <- update.logLik.Laplace(b, betas, sigma, gammas, gammas.bs, alpha, D)
    
    
    for (it in 1:iter) {
        # save parameter values in matrix
        Y.mat[it, ] <- c(betas, sigma)
        T.mat[it, ] <- c(gammas, gammas.bs, alpha)
        B.mat[it, ] <- D
        
        # linear predictors
        b.hat <- attr(new.b, "b")
        Zb <- rowSums(Z * b.hat[id, ], na.rm = TRUE)
        Zsb <- rowSums(Zs * b.hat[id.GK, ], na.rm = TRUE)
        eta.yx <- as.vector(X %*% betas)
        eta.yxT <- as.vector(Xtime %*% betas)
        Ys <- as.vector(Xs %*% betas) + Zsb
        eta.s <- as.vector(alpha * Ys)
        eta.ws <- as.vector(W2s %*% gammas.bs)

        sc1 <- - crossprod(X, y - eta.yx - Zb) / sigma^2
        sc2 <- - colSums(alpha * d * Xtime - exp(c(W1 %*% gammas)) * P * rowsum(exp(eta.ws + eta.s) * alpha * Xs, id.GK, reorder = FALSE))
        c(sc1 + sc2)
        
        
        
        
        # E-step -- compute EB estimates and traces
        for (i in 1:n) {
            # individual values
            eta.twi <- eta.tw[i]
            Ztime.i <- Ztime[i, ]
            # EB estimates, i.e., posterior modes
            bb <- attr(new.b, "b")[i, ]
            hes.b <- attr(new.b, "hes.b")[i, ]
            var.b <- attr(new.b, "vb")[i, ]            
            dim(var.b) <- dim(hes.b) <- c(ncz, ncz)
            # traces
            Ztime.b <- sum(Ztime.i * bb)
            eta.ti <- eta.twi + alpha * (eta.yxT[i] + Ztime.b)
            exp.eta.ti <- exp(eta.ti)
            trc1 <- - (alpha^3 * exp.eta.ti) * outer.Ztime[[i]]
            trc2 <- colSums(Ztime[i, ] * var.b)
            K <- var.b %*% trc1
            tr.var.b.trc1 <- - 0.5 * sum(diag(K))
            trc.y1[i, ] <- tr.var.b.trc1 * trc2
            trc.y2[i, ] <- - 0.5 * sum(- K * t(K)) * c(trc2 %o% trc2)
            L <- alpha * c(trc2 %o% trc2)
            M <- c(trc2 %o% colSums(Ztime.i * (K %*% var.b)))
            trc.y3[i, ] <- tr.var.b.trc1 * (L + M)
            ZtSZ <- c(crossprod(Ztime.i, solve(hes.b, Ztime.i)))
            P <- alpha * exp.eta.ti * ZtSZ - 1 / alpha
            trc.t1[i] <- tr.var.b.trc1 * P
            Q <- exp.eta.ti * ZtSZ * (alpha * Ztime.b + 1) - (alpha * Ztime.b + 2) / alpha^2
            trc.t2[i] <- tr.var.b.trc1 * Q
        }
        b <- attr(new.b, "b")
        hes.b <- attr(new.b, "hes.b")
        b.hat <- b + trc.y1
        vb.hat <- attr(new.b, "vb") + trc.y2 + trc.y3
        Zb <- rowSums(Z * b.hat[id, ])
        btZtZb <- drop(crossprod(Zb))
        outer.b.hat <- if (ncz == 1) apply(b.hat, 1, function (x) x %o% x) else t(apply(b.hat, 1, function (x) x %o% x))
        tr.tZZvarb <- sum(ZtZ * vb.hat)
        # compute log-likelihood and check convergence
        lgLik[it] <- as.vector(new.b)
        if (it > 2) {
            if (lgLik[it] < lgLik[it - 1]) {
                betas <- Y.mat[it - 1, 1:ncx]; sigma <- Y.mat[it - 1, ncx + 1]
                gammas <- T.mat[it - 1, 1:ncww]; alpha <- T.mat[it - 1, ncww + 1]
                D <- B.mat[it - 1,  ]
                if (!diag.D) dim(D) <- c(ncz, ncz)
                break
            } else {
                thets1 <- c(Y.mat[it - 1, ], T.mat[it - 1, ], B.mat[it - 1, ])
                thets2 <- c(Y.mat[it, ], T.mat[it, ], B.mat[it, ])
                check1 <- max(abs(thets2 - thets1) / (abs(thets1) + control$tol1)) < control$tol2
                check2 <- (lgLik[it] - lgLik[it - 1]) < control$tol3 * (abs(lgLik[it - 1]) + control$tol3)
                if (any(check1, check2)) {
                    conv <- TRUE
                    if (control$verbose) cat("\n\nconverged!\n")
                    break
                }
            }
        }
        
        # print results if verbose
        if (control$verbose) {
            cat("\n\niter:", it, "\n")
            cat("log-likelihood:", lgLik[it], "\n")
            cat("betas:", round(betas, 4), "\n")
            cat("sigma:", round(sigma, 4), "\n")
            cat("gammas:", round(gammas, 4), "\n")
            cat("alpha:", round(alpha, 4), "\n")
            cat("D:", if (!diag.D) round(D[lower.tri(D, TRUE)], 4) else round(D, 4), "\n")
        }
        
        # M-step
        mu <- y - eta.yx
        sigman <- sqrt(drop(crossprod(mu, mu - 2 * Zb) + btZtZb + tr.tZZvarb) / N)
        Dn <- matrix(colMeans(vb.hat + outer.b.hat), ncz, ncz)
        Dn <- if (diag.D) diag(Dn) else 0.5 * (Dn + t(Dn))
        Hbetas <- nearPD(fd.vec(betas, SclongCH))
        scbetas <- SclongCH(betas)
        betasn <- betas - c(solve(Hbetas, scbetas))
        thetas <- c(gammas, alpha)
        thetas[2:nk] <- log(diff(thetas[1:nk]))
        sc.thets <- ScsurvCH(thetas)
        # apply positive-definite modifications to the Hessian, if required
        H <- nearPD(cd.vec(thetas, ScsurvCH))
        step.thetas <- c(solve(H, sc.thets))
        # check if the step is too big, and scale if required
        sum.step <- sqrt(sum(step.thetas * step.thetas))
        sum.thet <- sqrt(sum(thetas * thetas))
        if (sum.step > (step.max <- control$step.max * max(sum.thet, ncww + 1)))
            step.thetas <- step.thetas * step.max / sum.step
        thetasn <- thetas - step.thetas
        # use backtracking if the log-likelihood has not increased
        new.b <- update.bCH(b, hes.b, betasn, sigman, thetasn, Dn, TRUE)
        dg0 <- - c(crossprod(sc.thets, step.thetas))
        g1 <- new.b
        g2 <- 0
        backtrack.step <- 0
        backtrack.fail <- FALSE
        while (!is.finite(g1) || g1 < lgLik[it] - 1e-04 * dg0) {
            if (!g2) {
                g2 <- g1
                lambda2 <- 1
                lambda1 <- 0.5 * dg0 / as.vector(g1 - lgLik[it] - dg0)
                if (is.na(lambda1) || (lambda1 > 1 | lambda1 < 0.1))
                    lambda1 <- 0.1
                g1 <- update.bCH(b, hes.b, betasn, sigman, thetas - lambda1 * step.thetas, Dn, TRUE)
            } else {
                L1 <- cbind(c(1 / lambda1^2, -lambda2 / lambda1^2), c(-1 / lambda2^2, lambda1 / lambda2^2))
                L2 <- c(g1 - dg0 * lambda1 - lgLik[it], g2 - dg0 * lambda2 - lgLik[it])
                k <- c(L1 %*% L2) / (lambda1 - lambda2)
                g2 <- g1
                lambda2 <- lambda1
                lambda1 <- (- k[2] + sqrt(k[2]^2 - 3 * k[1] * dg0)) / (3 * k[1])
                if (is.na(lambda1) || (lambda1 > 1 | lambda1 < 0.1))
                    lambda1 <- 0.1
                g1 <- update.bCH(b, hes.b, betasn, sigman, thetas - lambda1 * step.thetas, Dn, TRUE)
            }
            new.b <- g1
            thetasn <- thetas - lambda1 * step.thetas
            backtrack.step <- backtrack.step + 1
            if (backtrack.fail <- backtrack.step > control$backtrackSteps)
                break
        }
        if (control$verbose && backtrack.step > 0)
            cat("backtrack:", backtrack.step, "\tlambda:", lambda1, "\n")
        if (backtrack.fail) {
            if ((new.b <- update.bCH(b, hes.b, betas, sigma, thetasn, Dn, TRUE)) > lgLik[it] - 1e-04 * dg0) {
                betasn <- betas
                sigman <- sigma
            } else if ((new.b <- update.bCH(b, hes.b, betasn, sigman, thetasn, D, TRUE)) > lgLik[it] - 1e-04 * dg0) {
                Dn <- D
            } else {
                break
            }
        }
                
        # update parameter values
        betas <- betasn
        sigma <- sigman
        D <- Dn
        gammas <- thetasn[1:ncww]
        gammas[1:nk] <- cumsum(c(gammas[1], exp(gammas[2:nk])))
        alpha <- thetasn[ncww + 1]
    }
    thetsT <- c(gammas, alpha)
    thetsT[2:nk] <- log(diff(thetsT[1:nk]))
    thetas <- c(betas, log(sigma), thetsT, if (diag.D) log(D) else chol.transf(D))
    lgLik <- - LogLik.chLaplace(thetas, b = b)
    # if not converged, start quasi-Newton iterations
    if (!conv && !control$only.EM) {
        if (is.null(control$parscale))
            control$parscale <- rep(0.01, length(thetas))
        if (control$verbose)
            cat("\n\nquasi-Newton iterations start.\n\n")
        out <- if (control$optimizer == "optim") {
            optim(thetas, LogLik.chLaplace, Score.chLaplace, method = "BFGS",
                control = list(maxit = control$iter.qN, parscale = control$parscale, 
                trace = 10 * control$verbose), b = b)
        } else {
            nlminb(thetas, LogLik.chLaplace, Score.chLaplace, scale = control$parscale, 
                control = list(iter.max = control$iter.qN, trace = 1 * control$verbose), b = b)
        }
        if ((conv <- out$convergence) == 0 || - out[[2]] > lgLik) {
            lgLik <- - out[[2]]
            thetas <- out$par
            betas <- thetas[1:ncx]
            sigma <- exp(thetas[ncx + 1])
            gammas <- thetas[seq(ncx + 2, ncx + 1 + ncww)]
            gammas[1:nk] <- cumsum(c(gammas[1], exp(gammas[2:nk])))
            alpha <- thetas[ncx + ncww + 2]
            D <- thetas[seq(ncx + ncww + 3, length(thetas))]
            D <- if (diag.D) exp(D) else chol.transf(D)
            it <- it + if (control$optimizer == "optim") out$counts[1] else out$iterations
            # compute posterior moments for thetas after quasi-Newton
            new.b <- update.bCH(b, hes.b, betas, sigma, c(gammas, alpha), D)
            for (i in 1:n) {
                # individual values
                eta.twi <- eta.tw[i]
                Ztime.i <- Ztime[i, ]
                # EB estimates, i.e., posterior modes
                bb <- attr(new.b, "b")[i, ]
                hes.b <- attr(new.b, "hes.b")[i, ]
                var.b <- attr(new.b, "vb")[i, ]            
                dim(var.b) <- dim(hes.b) <- c(ncz, ncz)
                # traces
                Ztime.b <- sum(Ztime.i * bb)
                eta.ti <- eta.twi + alpha * (eta.yxT[i] + Ztime.b)
                exp.eta.ti <- exp(eta.ti)
                trc1 <- - (alpha^3 * exp.eta.ti) * outer.Ztime[[i]]
                trc2 <- colSums(Ztime[i, ] * var.b)
                K <- var.b %*% trc1
                tr.var.b.trc1 <- - 0.5 * sum(diag(K))
                trc.y1[i, ] <- tr.var.b.trc1 * trc2
                trc.y2[i, ] <- - 0.5 * sum(- K * t(K)) * c(trc2 %o% trc2)
                L <- alpha * c(trc2 %o% trc2)
                M <- c(trc2 %o% colSums(Ztime.i * (K %*% var.b)))
                trc.y3[i, ] <- tr.var.b.trc1 * (L + M)
                ZtSZ <- c(crossprod(Ztime.i, solve(hes.b, Ztime.i)))
                P <- alpha * exp.eta.ti * ZtSZ - 1 / alpha
                trc.t1[i] <- tr.var.b.trc1 * P
                Q <- exp.eta.ti * ZtSZ * (alpha * Ztime.b + 1) - (alpha * Ztime.b + 2) / alpha^2
                trc.t2[i] <- tr.var.b.trc1 * Q
            }
            b <- attr(new.b, "b")
            hes.b <- attr(new.b, "hes.b")
            b.hat <- b + trc.y1
            vb.hat <- attr(new.b, "vb") + trc.y2 + trc.y3
            Zb <- rowSums(Z * b.hat[id, ])
        }
    }
    # calculate Hessian matrix
    Hessian <- if (control$numeriDeriv == "fd") {
        fd.vec(thetas, Score.chLaplace, b = b, eps = control$eps.Hes)
    } else { 
        cd.vec(thetas, Score.chLaplace, b = b, eps = control$eps.Hes)
    }
    names(betas) <- names(initial.values$betas)
    if (!diag.D) dimnames(D) <- dimnames(initial.values$D) else names(D) <- names(initial.values$D)
    names(gammas) <- c(paste("bs.", 1:nk, sep = ""), colnames(W1))
    nams <- c(paste("Y.", c(names(betas), "sigma"), sep = ""), paste("T.", c(names(gammas), "alpha"), sep = ""),
        paste("B.", if (!diag.D) paste("D", seq(1, ncz * (ncz + 1) / 2), sep = "") else names(D), sep = "")
    )
    dimnames(Hessian) <- list(nams, nams)
    colnames(b.hat) <- colnames(x$Z)
    list(coefficients = list(betas = betas, sigma = sigma, gammas = gammas, alpha = alpha, D = as.matrix(D)), 
        Hessian = Hessian, logLik = lgLik, EB = list(post.b = b.hat, post.vb = vb.hat, Zb = Zb, 
        Ztimeb = rowSums(Ztime * b.hat)), knots = kn, iters = it, convergence = conv, n = n, N = N, ni = ni, d = d, 
        id = id)
}
