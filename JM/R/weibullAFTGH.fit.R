weibullAFTGH.fit <-
function (x, y, id, initial.values, scaleWB, parameterization, derivForm, control) {
    # response vectors
    logT <- as.vector(y$logT)
    d <- as.vector(y$d)
    y <- as.vector(y$y)
    # design matrices
    X <- x$X
    Xtime <- x$Xtime
    Xs <- x$Xs
    Xtime.deriv <- x$Xtime.deriv
    Xs.deriv <- x$Xs.deriv    
    Z <- x$Z
    Ztime <- x$Ztime
    Zs <- x$Zs
    Ztime.deriv <- x$Ztime.deriv
    Zs.deriv <- x$Zs.deriv
    W1 <- x$W
    WW <- if (is.null(W1)) as.matrix(rep(1, length(logT))) else cbind(1, W1)
    WintF.vl <- x$WintF.vl
    WintF.sl <- x$WintF.sl
    Ws.intF.vl <- x$Ws.intF.vl
    Ws.intF.sl <- x$Ws.intF.sl
    X <- dropAttr(X); Z <- dropAttr(Z); WW <- dropAttr(WW)
    WintF.vl <- dropAttr(WintF.vl); WintF.sl <- dropAttr(WintF.sl)
    Ws.intF.vl <- dropAttr(Ws.intF.vl); Ws.intF.sl <- dropAttr(Ws.intF.sl)
    if (parameterization == "value") {
        Xtime <- dropAttr(Xtime); Ztime <- dropAttr(Ztime); Xs <- dropAttr(Xs); Zs <- dropAttr(Zs)
    } else if (parameterization == "slope") {
        Xtime.deriv <- dropAttr(Xtime.deriv); Ztime.deriv <- dropAttr(Ztime.deriv)
        Xs.deriv <- dropAttr(Xs.deriv); Zs.deriv <- dropAttr(Zs.deriv)
    } else {
        Xtime <- dropAttr(Xtime); Ztime <- dropAttr(Ztime); Xs <- dropAttr(Xs); Zs <- dropAttr(Zs)
        Xtime.deriv <- dropAttr(Xtime.deriv); Ztime.deriv <- dropAttr(Ztime.deriv)
        Xs.deriv <- dropAttr(Xs.deriv); Zs.deriv <- dropAttr(Zs.deriv)
    }
    indFixed <- derivForm$indFixed
    indRandom <- derivForm$indRandom
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
    b <- as.matrix(expand.grid(rep(list(GH$x), ncz)))
    k <- nrow(b)
    wGH <- as.matrix(expand.grid(rep(list(GH$w), ncz)))    
    wGH <- 2^(ncz/2) * apply(wGH, 1, prod) * exp(rowSums(b * b))
    if (control$typeGH == "simple") {
        b <- sqrt(2) * t(control$inv.chol.VC %*% t(b))
        wGH <- wGH * control$det.inv.chol.VC
    } else { 
        b <- sqrt(2) * b
        VCdets <- control$det.inv.chol.VCs
    }
    dimnames(b) <- NULL
    b2 <- if (ncz == 1) b * b else t(apply(b, 1, function (x) x %o% x))
    Ztb <- Z %*% t(b)
    if (parameterization %in% c("value", "both")) {
        Ztime.b <- Ztime %*% t(b)
        Zsb <- Zs %*% t(b)
    }
    if (parameterization %in% c("slope", "both")) {
        if (length(indRandom) > 1 || indRandom) {
            Ztime.b.deriv <- Ztime.deriv %*% t(b[, indRandom, drop = FALSE])
            Zsb.deriv <- Zs.deriv %*% t(b[, indRandom, drop = FALSE])
        } else {
            Ztime.b.deriv <- matrix(0, nrow(Ztime.deriv), k)
            Zsb.deriv <- matrix(0, nrow(Zs.deriv), k)
        }
    }
    # Gauss-Kronrod rule
    st <- x$st
    log.st <- log(st)
    wk <- rep(x$wk, length(logT))
    P <- as.vector(x$P)
    id.GK <- rep(seq_along(logT), each = control$GKk)
    # pseudo-adaptive Gauss-Hermite
    if (control$typeGH != "simple") {
        lis.b <- vector("list", n)
        for (i in 1:n)
            lis.b[[i]] <- t(control$inv.chol.VCs[[i]] %*% t(b)) + 
                rep(control$ranef[i, ], each = k)
        lis.b2 <- lapply(lis.b, function (b) if (ncz == 1) b * b else
            t(apply(b, 1, function (x) x %o% x)))
        for (i in 1:n) {
            Ztb[id == i, ] <- Z[id == i, , drop = FALSE] %*% t(lis.b[[i]])
            if (parameterization %in% c("value", "both")) {
                Ztime.b[i, ] <- Ztime[i, , drop = FALSE] %*% t(lis.b[[i]])
                Zsb[id.GK == i, ] <- Zs[id.GK == i, ] %*% t(lis.b[[i]])
            }
            if (parameterization %in% c("slope", "both") && 
                    (length(indRandom) > 1 || indRandom)) {
                Ztime.b.deriv[i, ] <- Ztime.deriv[i, , drop = FALSE] %*% t(lis.b[[i]][, indRandom, drop = FALSE])
                Zsb.deriv[id.GK == i, ] <- Zs.deriv[id.GK == i, ] %*% t(lis.b[[i]][, indRandom, drop = FALSE])
            }            
        }
    }
    # initial values
    betas <- as.vector(initial.values$betas)
    sigma <- initial.values$sigma
    gammas <- as.vector(initial.values$gammas)
    alpha <- as.vector(initial.values$alpha)
    Dalpha <- as.vector(initial.values$Dalpha)
    sigma.t <- if (is.null(scaleWB)) initial.values$sigma.t else scaleWB
    D <- initial.values$D
    diag.D <- !is.matrix(D)
    if (!diag.D) dimnames(D) <- NULL else names(D) <- NULL
    # fix environments for functions
    environment(opt.survAFTWB) <- environment(gr.survAFTWB) <- environment()
    environment(opt.longAFTWB) <- environment(gr.longAFTWB) <- environment()
    environment(LogLik.weibullAFTGH) <- environment(Score.weibullAFTGH) <- environment()
    old <- options(warn = (-1))
    on.exit(options(old))
    # EM iterations
    iter <- control$iter.EM
    Y.mat <- matrix(0, iter + 1, ncx + 1)
    T.mat <- matrix(0, iter + 1, switch(parameterization, 
        "value" = ncww + 1 + ncol(WintF.vl),
        "slope" = ncww + 1 + ncol(WintF.sl),
        "both" = ncww + 1 + ncol(WintF.vl) + ncol(WintF.sl)))
    B.mat <- if (diag.D) matrix(0, iter + 1, ncz) else matrix(0, iter + 1, ncz * ncz)
    lgLik <- numeric(iter + 1)
    conv <- TRUE
    for (it in 1:iter) {
        # save parameter values in matrix
        Y.mat[it, ] <- c(betas, sigma)
        T.mat[it, ] <- switch(parameterization, "value" = c(gammas, alpha, sigma.t), 
            "slope" = c(gammas, Dalpha, sigma.t), "both" = c(gammas, alpha, Dalpha, sigma.t))
        B.mat[it,] <- D
        
        # linear predictors
        eta.yx <- as.vector(X %*% betas)
        eta.tw <- as.vector(WW %*% gammas)
        if (parameterization %in% c("value", "both")) {
            Y <- as.vector(Xtime %*% betas) + Ztime.b
            Ys <- as.vector(Xs %*% betas) + Zsb
            eta.t <- eta.tw + c(WintF.vl %*% alpha) * Y
            eta.s <- c(Ws.intF.vl %*% alpha) * Ys
        }
        if (parameterization %in% c("slope", "both")) {
            Y.deriv <- as.vector(Xtime.deriv %*% betas[indFixed]) + Ztime.b.deriv
            Ys.deriv <- as.vector(Xs.deriv %*% betas[indFixed]) + Zsb.deriv
            eta.t <- if (parameterization == "both")
                eta.t + c(WintF.sl %*% Dalpha) * Y.deriv
            else
                eta.tw + c(WintF.sl %*% Dalpha) * Y.deriv
            eta.s <- if (parameterization == "both")
                eta.s + c(Ws.intF.sl %*% Dalpha) * Ys.deriv
            else
                c(Ws.intF.sl %*% Dalpha) * Ys.deriv
        }
        
        # E-step
        mu.y <- eta.yx + Ztb
        logNorm <- dnorm(y, mu.y, sigma, TRUE)
        log.p.yb <- rowsum(logNorm, id, reorder = FALSE); dimnames(log.p.yb) <- NULL
        Vi <- exp(eta.tw) * P * rowsum(wk * exp(eta.s), id.GK, reorder = FALSE); dimnames(Vi) <- NULL
        log.hazard <- log(sigma.t) + (sigma.t - 1) * log(Vi) + eta.t
        log.survival <- - Vi^sigma.t
        log.p.tb <- d * log.hazard + log.survival
        log.p.b <- if (control$typeGH == "simple") {
            rep(dmvnorm(b, rep(0, ncz), D, TRUE), each = n)
        } else {
            matrix(dmvnorm(do.call(rbind, lis.b), rep(0, ncz), D, TRUE), n, k, byrow = TRUE)
        }
        p.ytb <- exp(log.p.yb + log.p.tb + log.p.b)
        if (control$typeGH != "simple")
            p.ytb <- p.ytb * VCdets
        p.yt <- c(p.ytb %*% wGH)
        p.byt <- p.ytb / p.yt
        post.b <- if (control$typeGH == "simple") {
            p.byt %*% (b * wGH)
        } else {
            sapply(seq_len(ncz), function (i)
                (p.byt * t(sapply(lis.b, "[", seq_len(k), i))) %*% wGH)
        }
        post.vb <- if (control$typeGH == "simple") { 
            if (ncz == 1) {
                c(p.byt %*% (b2 * wGH)) - c(post.b * post.b)
            } else {
                (p.byt %*% (b2 * wGH)) - t(apply(post.b, 1, function (x) x %o% x))
            }
        } else {
            dd <- sapply(seq_len(ncz^2), function (i)
                (p.byt * t(sapply(lis.b2, "[", seq_len(k), i))) %*% wGH)
            bb <- apply(post.b, 1, function (x) x %o% x)
            dd - if (ncz == 1) c(bb) else t(bb)
        }
        
        # compute log-likelihood
        log.p.yt <- log(p.yt)
        lgLik[it] <- sum(log.p.yt[is.finite(log.p.yt)])
        
        # print results if verbose
        if (control$verbose) {
            cat("\n\niter:", it, "\n")
            cat("log-likelihood:", lgLik[it], "\n")
            cat("betas:", round(betas, 4), "\n")
            cat("sigma:", round(sigma, 4), "\n")
            cat("gammas:", -round(gammas, 4), "\n")
            if (parameterization %in% c("value", "both"))
                cat("alpha:", -round(alpha, 4), "\n")
            if (parameterization %in% c("slope", "both"))
                cat("alphaD:", -round(Dalpha, 4), "\n")
            cat("sigma.t:", round(sigma.t, 4), "\n")
            cat("D:", if (!diag.D) round(D[lower.tri(D, TRUE)], 4) else round(D, 4), "\n")
        }
        
        # check convergence
        if (it > 5 && lgLik[it] > lgLik[it - 1]) {
            thets1 <- c(Y.mat[it - 1, ], T.mat[it - 1, ], B.mat[it - 1, ])
            thets2 <- c(Y.mat[it, ], T.mat[it, ], B.mat[it, ])
            check1 <- max(abs(thets2 - thets1) / (abs(thets1) + control$tol1)) < control$tol2
            check2 <- (lgLik[it] - lgLik[it - 1]) < control$tol3 * (abs(lgLik[it - 1]) + control$tol3)
            if (check1 || check2) {
                conv <- FALSE
                if (control$verbose)
                    cat("\n\nconverged!\ncalculating Hessian...\n")
                break
            }
        }
        if (iter == 0) break
        
        # M-step
        Zb <- rowSums(Z * post.b[id, ], na.rm = TRUE)
        mu <- y - eta.yx
        tr.tZZvarb <- sum(ZtZ * post.vb, na.rm = TRUE)        
        sigman <- sqrt(c(crossprod(mu, mu - 2 * Zb) + crossprod(Zb) + tr.tZZvarb) / N)
        Dn <- if (control$typeGH == "simple") {
            matrix(colMeans(p.byt %*% (b2 * wGH), na.rm = TRUE), ncz, ncz)
        } else {
            matrix(colMeans(dd, na.rm = TRUE), ncz, ncz)
        }
        Dn <- if (diag.D) diag(Dn) else 0.5 * (Dn + t(Dn))    
        Hbetas <- nearPD(fd.vec(betas, gr.longAFTWB))
        scbetas <- gr.longAFTWB(betas)
        betasn <- betas - c(solve(Hbetas, scbetas))
        list.thetas <- list(gammas = gammas, alpha = alpha, Dalpha = Dalpha, log.sigma.t = log(sigma.t))
        if (!is.null(scaleWB))
            list.thetas$log.sigma.t <- NULL
        list.thetas <- list.thetas[!sapply(list.thetas, is.null)]
        thetas <- unlist(as.relistable(list.thetas))
        optz.surv <- optim(thetas, opt.survAFTWB, gr.survAFTWB, method = "BFGS", 
            control = list(maxit = if (it < 5) 20 else 5, 
                parscale = if (it < 10) rep(0.01, length(thetas)) else rep(0.1, length(thetas))))
        thetasn <- relist(optz.surv$par, skeleton = list.thetas)

        # update parameter values
        betas <- betasn
        sigma <- sigman
        D <- Dn
        gammas <- thetasn$gammas
        alpha <- thetasn$alpha
        Dalpha <- thetasn$Dalpha
        sigma.t <- if (is.null(scaleWB)) exp(thetasn$log.sigma.t) else scaleWB        
    }
    list.thetas <- list(betas = betas, log.sigma = log(sigma), gammas = gammas, alpha = alpha, Dalpha = Dalpha,
        log.sigma.t = log(sigma.t), D = if (diag.D) log(D) else chol.transf(D))    
    if (!is.null(scaleWB))
        list.thetas$log.sigma.t <- NULL
    list.thetas <- list.thetas[!sapply(list.thetas, is.null)]
    thetas <- unlist(as.relistable(list.thetas))
    lgLik <- - LogLik.weibullAFTGH(thetas)
    
    # if not converged, start quasi-Newton iterations
    if (conv && !control$only.EM) {
        if (is.null(control$parscale))
            control$parscale <- rep(0.01, length(thetas))
        if (control$verbose)
            cat("\n\nquasi-Newton iterations start.\n\n")
        out <- if (control$optimizer == "optim") {
            optim(thetas, LogLik.weibullAFTGH, Score.weibullAFTGH, method = "BFGS",
                control = list(maxit = control$iter.qN, parscale = control$parscale, 
                trace = 10 * control$verbose))
        } else {
            nlminb(thetas, LogLik.weibullAFTGH, Score.weibullAFTGH, scale = control$parscale, 
                control = list(iter.max = control$iter.qN, trace = 1 * control$verbose))
        }
        if ((conv <- out$convergence) == 0 || - out[[2]] > lgLik) {
            lgLik <- - out[[2]]
            thetas <- relist(out$par, skeleton = list.thetas)
            betas <- thetas$betas
            sigma <- exp(thetas$log.sigma)
            gammas <- thetas$gammas
            alpha <- thetas$alpha
            Dalpha <- thetas$Dalpha
            sigma.t <- if (is.null(scaleWB)) exp(thetas$log.sigma.t) else scaleWB
            D <- thetas$D
            D <- if (diag.D) exp(D) else chol.transf(D)
            it <- it + if (control$optimizer == "optim") out$counts[1] else out$iterations
            # compute posterior moments for thetas after quasi-Newton
            eta.yx <- as.vector(X %*% betas)
            eta.tw <- as.vector(WW %*% gammas)
            if (parameterization %in% c("value", "both")) {
                Y <- as.vector(Xtime %*% betas) + Ztime.b
                Ys <- as.vector(Xs %*% betas) + Zsb
                eta.t <- eta.tw + c(WintF.vl %*% alpha) * Y
                eta.s <- c(Ws.intF.vl %*% alpha) * Ys
            }
            if (parameterization %in% c("slope", "both")) {
                Y.deriv <- as.vector(Xtime.deriv %*% betas[indFixed]) + Ztime.b.deriv
                Ys.deriv <- as.vector(Xs.deriv %*% betas[indFixed]) + Zsb.deriv
                eta.t <- if (parameterization == "both")
                    eta.t + c(WintF.sl %*% Dalpha) * Y.deriv
                else
                    eta.tw + c(WintF.sl %*% Dalpha) * Y.deriv
                eta.s <- if (parameterization == "both")
                    eta.s + c(Ws.intF.sl %*% Dalpha) * Ys.deriv
                else
                    c(Ws.intF.sl %*% Dalpha) * Ys.deriv
            }
            exp.eta.tw <- exp(eta.tw)
            mu.y <- eta.yx + Ztb
            logNorm <- dnorm(y, mu.y, sigma, TRUE)
            log.p.yb <- rowsum(logNorm, id)
            Vi <- exp(eta.tw) * P * rowsum(wk * exp(eta.s), id.GK, reorder = FALSE); dimnames(Vi) <- NULL
            log.hazard <- log(sigma.t) + (sigma.t - 1) * log(Vi) + eta.t
            log.survival <- - Vi^sigma.t
            log.p.tb <- d * log.hazard + log.survival            
            log.p.b <- if (control$typeGH == "simple") {
                rep(dmvnorm(b, rep(0, ncz), D, TRUE), each = n)
            } else {
                matrix(dmvnorm(do.call(rbind, lis.b), rep(0, ncz), D, TRUE), n, k, byrow = TRUE)
            }
            p.ytb <- exp(log.p.yb + log.p.tb + log.p.b)
            if (control$typeGH != "simple")
                p.ytb <- p.ytb * VCdets
            p.yt <- c(p.ytb %*% wGH)
            p.byt <- p.ytb / p.yt
            post.b <- if (control$typeGH == "simple") {
                p.byt %*% (b * wGH)
            } else {
                sapply(seq_len(ncz), function (i)
                    (p.byt * t(sapply(lis.b, "[", seq_len(k), i))) %*% wGH)
            }
            post.vb <- if (control$typeGH == "simple") { 
                if (ncz == 1) {
                    c(p.byt %*% (b2 * wGH)) - c(post.b * post.b)
                } else {
                    (p.byt %*% (b2 * wGH)) - t(apply(post.b, 1, function (x) x %o% x))
                }
            } else {
                dd <- sapply(seq_len(ncz^2), function (i)
                    (p.byt * t(sapply(lis.b2, "[", seq_len(k), i))) %*% wGH)
                bb <- apply(post.b, 1, function (x) x %o% x)
                dd - if (ncz == 1) c(bb) else t(bb)
            }
            Zb <- if (ncz == 1) post.b[id] else rowSums(Z * post.b[id, ], na.rm = TRUE)
        }
    }
    # calculate Score vector
    Score <- Score.weibullAFTGH(unlist(thetas))
    # calculate Hessian matrix
    Hessian <- if (control$numeriDeriv == "fd") {
        fd.vec(unlist(thetas), Score.weibullAFTGH, eps = control$eps.Hes)
    } else { 
        cd.vec(unlist(thetas), Score.weibullAFTGH, eps = control$eps.Hes)
    }    
    names(betas) <- names(initial.values$betas)
    if (!diag.D) dimnames(D) <- dimnames(initial.values$D) else names(D) <- names(initial.values$D)
    names(gammas) <- c("(Intercept)", colnames(W1))
    nm.alph <- colnames(x$WintF.vl)
    nm.alph <- if (!is.null(nm.alph)) {
        if (nm.alph[1] == "(Intercept)")
            c("", nm.alph[-1])
        else
            nm.alph
    } else {
        "alpha"
    }
    nm.Dalph <- colnames(x$WintF.sl)
    nm.Dalph <- if (!is.null(nm.Dalph)) {
        if (nm.Dalph[1] == "(Intercept)")
            c("", nm.Dalph[-1])
        else
            nm.Dalph
    } else {
        "alpha.s"
    }
    gg <- switch(parameterization, "value" = nm.alph, "slope" = nm.Dalph, "both" = c(nm.alph, nm.Dalph))
    if (parameterization %in% c("value", "both"))
        names(alpha) <- nm.alph
    if (parameterization %in% c("slope", "both"))
        names(Dalpha) <- nm.Dalph
    nams <- c(paste("Y.", c(names(betas), "sigma"), sep = ""), 
        paste("T.", if (is.null(scaleWB)) c(names(gammas), gg, "sigma.t") else c(names(gammas), gg), sep = ""),
        paste("B.", if (!diag.D) paste("D", seq(1, ncz * (ncz + 1) / 2), sep = "") else names(D), sep = ""))
    dimnames(Hessian) <- list(nams, nams)
    colnames(post.b) <- colnames(x$Z)
    list(coefficients = list(betas = betas, sigma = sigma, gammas = gammas, alpha = alpha, Dalpha = Dalpha, sigma.t = sigma.t, 
        D = as.matrix(D)), Score = Score, Hessian = Hessian, logLik = lgLik, EB = list(post.b = post.b, post.vb = post.vb, 
        Zb = if (iter == 0) rowSums(Z * post.b[id, ], na.rm = TRUE) else Zb, 
        Ztimeb = if (parameterization %in% c("value", "both")) rowSums(Ztime * post.b) else NULL,
        Ztimeb.deriv = if (parameterization %in% c("slope", "both")) {
            if (indRandom) rowSums(Ztime.deriv * post.b[, indRandom, drop = FALSE]) else rep(0, nrow(Ztime.deriv))
        } else NULL), 
        iters = it, convergence = conv, n = n, N = N, ni = ni, d = d, id = id, scaleWB = scaleWB)
}
