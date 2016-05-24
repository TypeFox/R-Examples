piecewiseAFTGH.fit <-
function (x, y, id, initial.values, control) {
    # response vectors
    logT <- as.vector(y$logT)
    d <- as.vector(y$d)
    y <- as.vector(y$y)
    # design matrices
    X <- x$X
    Xtime <- x$Xtime
    Xs <- x$Xs
    Z <- x$Z
    Ztime <- x$Ztime
    Zs <- x$Zs
    W1 <- x$W
    WW <- if (is.null(W1)) as.matrix(rep(1, length(logT))) else cbind(1, W1)
    dimnames(X) <- dimnames(Xtime) <- dimnames(Xs) <- dimnames(Z) <- dimnames(Ztime) <- dimnames(Zs) <- dimnames(WW) <- NULL
    attr(X, "assign") <- attr(X, "contrasts") <- attr(Xtime, "assign") <- attr(Xtime, "contrasts") <- NULL
    attr(Xs, "assign") <- attr(Xs, "contrasts") <- attr(Zs, "assign") <- attr(Zs, "contrasts") <- attr(Z, "assign") <- attr(Ztime, "assign") <- NULL
    # sample size settings
    ncx <- ncol(X)
    ncz <- ncol(Z)
    ncww <- ncol(WW)
    n <- length(logT)
    N <- length(y)
    ni <- as.vector(tapply(id, id, length))
    # crossproducts and others
    XtX <- crossprod(X)
    ZtZ <- lapply(split(Z, id), function (x) crossprod(matrix(x, ncol = ncz)))
    names(ZtZ) <- NULL
    ZtZ <- matrix(unlist(ZtZ), n, ncz * ncz, TRUE)
    outer.Ztime <- lapply(1:n, function (x) Ztime[x, ] %o% Ztime[x, ])
    # Gauss-Hermite quadrature rule components
    GH <- gauher(control$GHk)
    b <- as.matrix(expand.grid(lapply(1:ncz, function (k, u) u$x, u = GH)))
    k <- nrow(b)
    wGH <- as.matrix(expand.grid(lapply(1:ncz, function (k, u) u$w, u = GH)))
    wGH <- 2^(ncz/2) * apply(wGH, 1, prod) * exp(rowSums(b * b)) * control$det.inv.chol.VC
    b <- sqrt(2) * t(control$inv.chol.VC %*% t(b)); dimnames(b) <- NULL
    b2 <- if (ncz == 1) b * b else t(apply(b, 1, function (x) x %o% x))
    Ztb <- Z %*% t(b)
    Ztime.b <- Ztime %*% t(b)
    Zsb <- Zs %*% t(b)
    # Gauss-Kronrod rule
    st <- x$st
    log.st <- log(st)
    wk <- rep(x$wk, length(logT))
    P <- as.vector(x$P)
    id.GK <- rep(seq_along(logT), each = nk)
    # initial values
    betas <- as.vector(initial.values$betas)
    sigma <- initial.values$sigma
    gammas <- as.vector(initial.values$gammas)
    alpha <- as.vector(initial.values$alpha)
    sigma.t <- initial.values$sigma.t
    D <- initial.values$D
    diag.D <- !is.matrix(D)
    if (!diag.D) dimnames(D) <- NULL else names(D) <- NULL
    # fix environments for functions
    environment(opt.survAFTPC) <- environment(gr.survAFTPC) <- environment()
    environment(opt.longAFTPC) <- environment(gr.longAFTPC) <- environment()
    environment(LogLik.piecewiseAFTGH) <- environment(Score.piecewiseAFTGH) <- environment()
    old <- options(warn = (-1))
    on.exit(options(old))
    # EM iterations
    iter <- control$iter.EM
    Y.mat <- matrix(0, iter, ncx + 1)
    T.mat <- matrix(0, iter, ncww + 2)
    B.mat <- if (diag.D) matrix(0, iter, ncz) else matrix(0, iter, ncz * ncz)
    lgLik <- numeric(iter)
    conv <- FALSE
    for (it in 1:iter) {
        # save parameter values in matrix
        Y.mat[it, ] <- c(betas, sigma)
        T.mat[it, ] <- c(gammas, alpha, sigma.t)
        B.mat[it,] <- D
        
        # linear predictors
        eta.yx <- as.vector(X %*% betas)
        eta.yxT <- as.vector(Xtime %*% betas)
        eta.tw <- as.vector(WW %*% gammas)
        Y <- eta.yxT + Ztime.b
        Ys <- as.vector(Xs %*% betas) + Zsb
        eta.t <- eta.tw + alpha * Y
        eta.s <- alpha * Ys
        
        # E-step
        mu.y <- eta.yx + Ztb
        logNorm <- dnorm(y, mu.y, sigma, TRUE)
        log.p.yb <- rowsum(logNorm, id, reorder = FALSE); dimnames(log.p.yb) <- NULL
        Vi <- exp(eta.tw) * P * rowsum(wk * exp(eta.s), id.GK, reorder = FALSE); dimnames(Vi) <- NULL
        log.hazard <- log(sigma.t) + (sigma.t - 1) * log(Vi) + eta.t
        log.survival <- - Vi^sigma.t
        log.p.tb <- d * log.hazard + log.survival
        log.p.b <- if (ncz == 1) {
            dnorm(b, sd = sqrt(D), log = TRUE)
        } else {
            if (diag.D) {
                rowSums(dnorm(b, sd = rep(sqrt(D), each = k), log = TRUE))
            } else {
                dmvnorm(b, rep(0, ncz), D, TRUE)
            }
        }
        p.ytb <- exp((log.p.yb + log.p.tb) + rep(log.p.b, each = n))
        p.yt <- c(p.ytb %*% wGH)
        p.byt <- p.ytb / p.yt
        post.b <- p.byt %*% (b * wGH)
        post.vb <- if (ncz == 1) {
            c(p.byt %*% (b2 * wGH)) - c(post.b * post.b)
        } else {
            (p.byt %*% (b2 * wGH)) - t(apply(post.b, 1, function (x) x %o% x))
        }
        
        # compute log-likelihood
        log.p.yt <- log(p.yt)
        lgLik[it] <- sum(log.p.yt[is.finite(log.p.yt)], na.rm = TRUE)
        
        # print results if verbose
        if (control$verbose) {
            cat("\n\niter:", it, "\n")
            cat("log-likelihood:", lgLik[it], "\n")
            cat("betas:", round(betas, 4), "\n")
            cat("sigma:", round(sigma, 4), "\n")
            cat("gammas:", -round(gammas, 4), "\n")
            cat("alpha:", -round(alpha, 4), "\n")
            cat("sigma.t:", round(sigma.t, 4), "\n")
            cat("D:", if (!diag.D) round(D[lower.tri(D, TRUE)], 4) else round(D, 4), "\n")
        }
        
        # check convergence
        if (it > 5) {
            if (lgLik[it] < lgLik[it - 1]) {
                betas <- Y.mat[it - 1, 1:ncx]
                sigma <- Y.mat[it - 1, ncx + 1]
                gammas <- T.mat[it - 1, 1:ncww]
                alpha <- T.mat[it - 1, ncww + 1]
                sigma.t <- T.mat[it - 1, ncww + 2]
                D <- B.mat[it - 1,  ]
                if (!diag.D) dim(D) <- c(ncz, ncz)
                break
            } else {
                thets1 <- c(Y.mat[it - 1, ], T.mat[it - 1, ], B.mat[it - 1, ])
                thets2 <- c(Y.mat[it, ], T.mat[it, ], B.mat[it, ])
                check1 <- max(abs(thets2 - thets1) / (abs(thets1) + control$tol1)) < control$tol2
                check2 <- (lgLik[it] - lgLik[it - 1]) < control$tol3 * (abs(lgLik[it - 1]) + control$tol3)
                if (check1 || check2) {
                    conv <- TRUE
                    if (control$verbose) cat("\n\nconverged!\ncalculating Hessian...\n")
                    break
                }
            }
        }
                
        # M-step
        Zb <- rowSums(Z * post.b[id, ], na.rm = TRUE)
        mu <- y - eta.yx
        tr.tZZvarb <- sum(ZtZ * post.vb, na.rm = TRUE)
        sigman <- sqrt(c(crossprod(mu, mu - 2 * Zb) + crossprod(Zb) + tr.tZZvarb) / N)
        Dn <- matrix(colMeans(p.byt %*% (b2 * wGH), na.rm = TRUE), ncz, ncz)
        Dn <- if (diag.D) diag(Dn) else 0.5 * (Dn + t(Dn))
        Hbetas <- nearPD(fd.vec(betas, gr.longAFTPC))
        scbetas <- gr.longAFTPC(betas)
        betasn <- betas - c(solve(Hbetas, scbetas))
        thetas <- c(gammas, alpha, log(sigma.t))
        optz.surv <- optim(thetas, opt.survAFTPC, gr.survAFTPC, method = "BFGS", 
            control = list(maxit = if (it < 5) 20 else 5, 
                parscale = if (it < 10) rep(0.01, length(thetas)) else rep(0.1, length(thetas))))
        thetasn <- optz.surv$par

        # update parameter values
        betas <- betasn
        sigma <- sigman
        D <- Dn
        gammas <- thetasn[1:ncww]
        alpha <- thetasn[ncww + 1]
        sigma.t <- exp(thetasn[ncww + 2])
    }
    thetas <- c(betas, log(sigma), gammas, alpha, log(sigma.t), if (diag.D) log(D) else chol.transf(D))
    lgLik <- - LogLik.piecewiseAFTGH(thetas)    
    # if not converged, start quasi-Newton iterations
    if (!conv && !control$only.EM) {
        if (is.null(control$parscale))
            control$parscale <- rep(0.01, length(thetas))
        if (control$verbose)
            cat("\n\nquasi-Newton iterations start.\n\n")
        out <- if (control$optimizer == "optim") {
            optim(thetas, LogLik.piecewiseAFTGH, Score.piecewiseAFTGH, method = "BFGS",
                control = list(maxit = control$iter.qN, parscale = control$parscale, 
                trace = 10 * control$verbose))
        } else {
            nlminb(thetas, LogLik.piecewiseAFTGH, Score.piecewiseAFTGH, scale = control$parscale, 
                control = list(iter.max = control$iter.qN, trace = 1 * control$verbose))
        }
        if ((conv <- out$convergence) == 0 || - out[[2]] > lgLik) {
            lgLik <- - out[[2]]
            thetas <- out$par
            betas <- thetas[1:ncx]
            sigma <- exp(thetas[ncx + 1])
            gammas <- thetas[seq(ncx + 2, ncx + 1 + ncww)]
            alpha <- thetas[ncx + ncww + 2]
            sigma.t <- exp(thetas[ncx + ncww + 3])
            D <- thetas[seq(ncx + ncww + 4, length(thetas))]
            D <- if (diag.D) exp(D) else chol.transf(D)
            it <- it + if (control$optimizer == "optim") out$counts[1] else out$iterations
            # compute posterior moments for thetas after quasi-Newton
            eta.yx <- as.vector(X %*% betas)
            eta.yxT <- as.vector(Xtime %*% betas)
            eta.tw <- as.vector(WW %*% gammas)
            exp.eta.tw <- exp(eta.tw)
            Y <- eta.yxT + Ztime.b
            Ys <- as.vector(Xs %*% betas) + Zsb
            eta.t <- eta.tw + alpha * Y
            eta.s <- alpha * Ys
            mu.y <- eta.yx + Ztb
            logNorm <- dnorm(y, mu.y, sigma, TRUE)
            log.p.yb <- rowsum(logNorm, id)
            Vi <- exp(eta.tw) * P * rowsum(wk * exp(eta.s), id.GK, reorder = FALSE); dimnames(Vi) <- NULL
            log.hazard <- log(sigma.t) + (sigma.t - 1) * log(Vi) + eta.t
            log.survival <- - Vi^sigma.t
            log.p.tb <- d * log.hazard + log.survival            
            log.p.b <- if (ncz == 1) {
                dnorm(b, sd = sqrt(D), log = TRUE)
            } else {
                if (diag.D) {
                    rowSums(dnorm(b, sd = rep(sqrt(D), each = k), log = TRUE))
                } else {
                    dmvnorm(b, rep(0, ncz), D, TRUE)
                }
            }
            p.ytb <- exp((log.p.yb + log.p.tb) + rep(log.p.b, each = n))
            dimnames(p.ytb) <- NULL
            p.yt <- c(p.ytb %*% wGH)
            p.byt <- p.ytb / p.yt
            post.b <- p.byt %*% (b * wGH)
            post.vb <- if (ncz == 1) {
                c(p.byt %*% (b2 * wGH)) - c(post.b * post.b)
            } else {
                (p.byt %*% (b2 * wGH)) - t(apply(post.b, 1, function (x) x %o% x))
            }
            Zb <- if (ncz == 1) post.b[id] else rowSums(Z * post.b[id, ], na.rm = TRUE)
        }
    }
    # calculate Hessian matrix
    Hessian <- if (control$numeriDeriv == "fd") {
        fd.vec(thetas, Score.piecewiseAFTGH, eps = control$eps.Hes)
    } else { 
        cd.vec(thetas, Score.piecewiseAFTGH, eps = control$eps.Hes)
    }
    names(betas) <- names(initial.values$betas)
    if (!diag.D) dimnames(D) <- dimnames(initial.values$D) else names(D) <- names(initial.values$D)
    names(gammas) <- c("(Intercept)", colnames(W1))
    nams <- c(paste("Y.", c(names(betas), "sigma"), sep = ""), paste("T.", c(names(gammas), "alpha", "sigma.t"), sep = ""),
        paste("B.", if (!diag.D) paste("D", seq(1, ncz * (ncz + 1) / 2), sep = "") else names(D), sep = ""))
    dimnames(Hessian) <- list(nams, nams)
    colnames(post.b) <- colnames(x$Z)
    list(coefficients = list(betas = betas, sigma = sigma, gammas = gammas, alpha = alpha, sigma.t = sigma.t, 
        D = as.matrix(D)), Hessian = Hessian, logLik = lgLik, EB = list(post.b = post.b, post.vb = post.vb, Zb = Zb, 
        Ztimeb = rowSums(Ztime * post.b)), iters = it, convergence = conv, n = n, N = N, ni = ni, d = d, id = id)
}
