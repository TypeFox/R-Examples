rocJM <-
function (object, dt, data, idVar = "id", directionSmaller = NULL, cc = NULL, 
        min.cc = NULL, max.cc = NULL, optThr = c("sens*spec", "youden"), 
        diffType = c("absolute", "relative"), 
        abs.diff = 0, rel.diff = 1, M = 300, burn.in = 100, scale = 1.6) {
    if (!inherits(object, "jointModel"))
        stop("Use only with 'jointModel' objects.\n")
    if (object$CompRisk)
        stop("rocJM() is not currently implemented for ",
            "competing risks joint models.\n")
    if (!is.data.frame(data) || nrow(data) == 0)
        stop("'data' must be a data.frame with more than one rows.\n")
    if (is.null(data[[idVar]]))
        stop("'idVar' not in 'data.\n'")
    optThr  <- match.arg(optThr)
    method <- object$method
    if (is.null(directionSmaller)) {
        directionSmaller <- as.logical(if (method == "weibull-AFT-GH")
            object$coefficients$alpha > 0 else object$coefficients$alpha < 0)
    }
    if (!length(directionSmaller)) {
        stop("you need to define the 'directionSmaller' argument; ", 
            "see the online help file for more details.\n")
    }
    if (is.null(cc) || !is.numeric(cc)) {
        pc <- quantile(object$y$y, c(0.05, 0.95), names = FALSE)
        if (is.null(min.cc) || !is.numeric(min.cc))
            min.cc <- min(min.cc, min(object$y$y) - pc[1])
        if (is.null(max.cc) || !is.numeric(max.cc))
            max.cc <- max(max.cc, max(object$y$y) + pc[2])
        cc <- seq(min.cc, max.cc, length.out = 251)
    }
    diffType <- match.arg(diffType)
    if (diffType == "absolute") {
        lag <- length(abs.diff)
        cc <- matrix(cc, length(cc), lag) + 
            rep(abs.diff, each = length(cc))
    } else {
        lag <- length(rel.diff)
        vv <- outer(cc, rel.diff - 1)
        #if (directionSmaller) {
        #    outer(cc, rel.diff - 1)
        #} else {
        #    outer(cc, rel.diff - 1)
        #}
        cc <- matrix(cc, length(cc), lag) + sign(cc) * vv
        #cc <- matrix(cc, length(cc), lag) * rep(rel.diff, each = length(cc))
    }
    timeVar <- object$timeVar
    interFact <- object$interFact
    parameterization <- object$parameterization
    derivForm <- object$derivForm
    indFixed <- derivForm$indFixed
    indRandom <- derivForm$indRandom
    LongFormat <- object$LongFormat
    id <- as.numeric(unclass(data[[idVar]]))
    id <- match(id, unique(id))
    TermsX <- delete.response(object$termsYx)
    TermsZ <- object$termsYz
    TermsX.deriv <- object$termsYx.deriv
    TermsZ.deriv <- object$termsYz.deriv
    mfX <- model.frame(TermsX, data = data)
    mfZ <- model.frame(TermsZ, data = data)
    formYx <- reformulate(attr(TermsX, "term.labels"))
    formYz <- object$formYz
    X <- model.matrix(formYx, mfX)
    Z <- model.matrix(formYz, mfZ)
    TermsT <- object$termsT
    data.id <- if (LongFormat) data else data[!duplicated(id), ]
    data.id <- if (LongFormat) {
        nams.ind <- all.vars(delete.response(TermsT))
        ind <- !duplicated(data[nams.ind])
        data[ind, ]
    } else data[!duplicated(id), ]
    idT <- data.id[[idVar]]
    idT <- match(idT, unique(idT))
    mfT <- model.frame(delete.response(TermsT), data = data.id)
    formT <- if (!is.null(kk <- attr(TermsT, "specials")$strata)) {
        strt <- eval(attr(TermsT, "variables"), data.id)[[kk]]
        tt <- drop.terms(TermsT, kk - 1, keep.response = FALSE)
        reformulate(attr(tt, "term.labels"))
    } else if (!is.null(kk <- attr(TermsT, "specials")$cluster)) {
        tt <- drop.terms(TermsT, kk - 1, keep.response = FALSE)
        reformulate(attr(tt, "term.labels"))
    } else {
        tt <- attr(delete.response(TermsT), "term.labels")
        if (length(tt)) reformulate(tt) else reformulate("1")
    }
    W <- model.matrix(formT, mfT)
    WintF.vl <- WintF.sl <- as.matrix(rep(1, nrow(data.id)))
    if (!is.null(interFact)) {
        if (!is.null(interFact$value))
            WintF.vl <- model.matrix(interFact$value, data = data.id)
        if (!is.null(interFact$slope))
            WintF.sl <- model.matrix(interFact$slope, data = data.id)
    }
    obs.times <- split(data[[timeVar]], id)
    last.time <- tapply(data[[timeVar]], id, tail, n = 1)
    tDt <- lapply(last.time, function (t) {s <- t + dt; s[s > t]})
    n <- object$n
    n.tp <- length(last.time)
    ncx <- ncol(X)
    ncz <- ncol(Z)
    ncww <- ncol(W)
    betas <- object$coefficients$betas
    sigma <- object$coefficients$sigma
    D <- object$coefficients$D
    diag.D <- ncol(D) == 1 & nrow(D) > 1
    D <- if (diag.D) diag(c(D)) else D
    gammas <- object$coefficients$gammas
    alpha <- object$coefficients$alpha
    Dalpha <- object$coefficients$Dalpha
    sigma.t <- object$coefficients$sigma.t
    xi <- object$coefficients$xi
    gammas.bs <- object$coefficients$gammas.bs
    list.thetas <- list(betas = betas, log.sigma = log(sigma), 
        gammas = gammas, alpha = alpha, Dalpha = Dalpha, 
        log.sigma.t = if (is.null(sigma.t)) NULL else log(sigma.t), 
        log.xi = if (is.null(xi)) NULL else log(xi), gammas.bs = gammas.bs, 
        D = if (diag.D) log(diag(D)) else chol.transf(D))
    if (method %in% c("weibull-PH-GH", "weibull-AFT-GH") && 
            !is.null(object$scaleWB)) {
        list.thetas$log.sigma.t <- NULL
    }
    if (method %in% c("piecewise-PH-GH", "spline-PH-GH")) {
        if (ncww == 1) {
            W <- NULL
            ncww <- 0
        } else {
            W <- W[, -1, drop = FALSE]
            ncww <- ncww - 1
        }
        Q <- object$x$Q
    }
    list.thetas <- list.thetas[!sapply(list.thetas, is.null)]
    if (!method %in% c("weibull-PH-GH", "weibull-AFT-GH", 
            "piecewise-PH-GH", "spline-PH-GH")) {
        stop("\nrocJM() is not yet available for this type of joint model.")
    }
    list.thetas <- list.thetas[!sapply(list.thetas, is.null)]
    thetas <- unlist(as.relistable(list.thetas))
    Var.thetas <- vcov(object)
    environment(log.posterior.b) <- environment(S.b) <- environment()
    environment(ModelMats) <- environment()
    # Simulate
    out <- vector("list", M)
    success.rate <- matrix(FALSE, M, n.tp)
    b.old <- b.new <- mvrnorm(n.tp, rep(0, ncz), D)
    if (n.tp == 1)
        dim(b.old) <- dim(b.new) <- c(1, ncz)
    y.new <- numeric(nrow(data))
    survTimes <- c(outer(dt, last.time, "+"))
    # construct model matrices to calculate the survival functions
    survMats <- survMats.last <- vector("list", n.tp)
    for (i in seq_len(n.tp)) {
        survMats[[i]] <- lapply(tDt[[i]], ModelMats, ii = i,
            obs.times = obs.times, survTimes = survTimes)
        survMats.last[[i]] <- ModelMats(last.time[i], ii = i,
            obs.times = obs.times, survTimes = survTimes) 
    }
    for (m in 1:M) {
        # Step 1: simulate new parameter values
        thetas.new <- mvrnorm(1, thetas, Var.thetas)
        thetas.new <- relist(thetas.new, skeleton = list.thetas)
        betas.new <- thetas.new$betas
        sigma.new <- exp(thetas.new$log.sigma)
        gammas.new <- thetas.new$gammas
        alpha.new <- thetas.new$alpha
        Dalpha.new <- thetas.new$Dalpha
        D.new <- thetas.new$D
        D.new <- if (diag.D) exp(D.new) else chol.transf(D.new)
        if (method == "weibull-PH-GH" || method == "weibull-AFT-GH") {
            sigma.t.new <- if (is.null(object$scaleWB)) 
                exp(thetas.new$log.sigma.t)
            else
                object$scaleWB
        } else if (method == "piecewise-PH-GH") {
            xi.new <- exp(thetas.new$log.xi)
        } else if (method == "spline-PH-GH") {
            gammas.bs.new <- thetas.new$gammas.bs
        }
        pi.dt.t <- vector("list", n.tp)
        F.y.t <- matrix(0, n.tp, nrow(cc))
        for (i in seq_len(n.tp)) {
            # Step 2: simulate new y values
            id.i <- id %in% i
            x.i <- X[id.i, , drop = FALSE]
            z.i <- Z[id.i, , drop = FALSE]
            mu.new <- as.vector(x.i %*% betas.new) + rowSums(z.i * rep(b.new[i, ], each = nrow(z.i)))
            y.new[id.i] <- rnorm(length(mu.new), mu.new, sigma.new)
            # Step 3: simulate new random effects values
            if (length(tDt[[i]])) {
                ff <- function (b, y, tt, mm, i) -log.posterior.b(b, y = y, Mats = tt, method = mm, ii = i)
                gg <- function (b, y, tt, mm, i) cd(b, ff, y = y, tt = tt, mm = mm, i = i)
                opt <- optim(rep(0, ncz), ff, gg, y = y.new, tt = survMats.last, mm = method, i = i, 
                    method = "BFGS", hessian = TRUE)
                mode.b <- opt$par
                var.b <- scale * solve(opt$hessian)
                proposed.b <- rmvt(1, mode.b, var.b, 4)
                dmvt.old <- dmvt(b.old[i, ], mode.b, var.b, 4, TRUE)
                dmvt.proposed <- dmvt(proposed.b, mode.b, var.b, 4, TRUE)
                a <- min(exp(log.posterior.b(proposed.b, y.new, survMats.last, method, ii = i) + dmvt.old - 
                    log.posterior.b(b.old[i, ], y.new, survMats.last, method, ii = i) - dmvt.proposed), 1)
                ind <- runif(1) <= a
                success.rate[m, i] <- ind
                if (ind)
                    b.new[i, ] <- proposed.b
                # Step 4a: compute Pr(T > t_k | T > t_{k - 1}, b.new; theta.new)
                S.last <- S.b(last.time[i], b.new[i, ], i, survMats.last[[i]])
                S.pred <- numeric(length(tDt[[i]]))
                for (l in seq_along(S.pred))
                    S.pred[l] <- S.b(tDt[[i]][l], b.new[i, ], i, survMats[[i]][[l]])
                pi.dt.t[[i]] <- S.pred / S.last
                # Step 4b: compute Pr(y(t) <= c | b.new; theta.new)
                mu.new <- as.vector(x.i %*% betas.new) + 
                    rowSums(z.i * rep(b.new[i, ], each = nrow(z.i)))                
                ccc <- cc[, seq(1, min(nrow(z.i), ncol(cc))), drop = FALSE]
                ppp <- pnorm(ccc, mean = rep(tail(mu.new, lag), each = nrow(cc)), 
                    sd = sigma.new, log.p = TRUE, lower.tail = directionSmaller)
                F.y.t[i, ] <- exp(rowSums(ppp))
             } else {
                pi.dt.t[[i]] <- F.y.t[i, ] <- NA
            }
        }
        b.old <- b.new
        out[[m]] <- c(pi.dt.t = list(pi.dt.t), F.y.t = list(F.y.t))
    }
    res <- vector("list", n.tp)
    for (i in seq_len(n.tp)) {
        pp <- sapply(out, function (x) x[[1]][[i]])
        if (!is.matrix(pp))
            pp <- rbind(pp)
        pp <- pp[, -seq_len(burn.in), drop = FALSE]
        ff <- sapply(out, function (x) x[[2]][i, ])[, -seq_len(burn.in)]
        M <- ncol(ff)
        p1 <- rowMeans(pp, na.rm = TRUE)
        v1 <- apply(t(pp), 2, sd, na.rm = TRUE)^2 / M
        rr2 <- vector("list", M)
        for (m in seq_len(M))
            rr2[[m]] <- outer(1 - pp[, m], ff[, m])
        p2 <- Reduce("+", rr2) / M
        v2 <- apply(array(unlist(rr2), c(length(dt), nrow(cc), length(rr2))), c(1, 2), sd)^2 / M
        rr3 <- vector("list", M)
        for (m in seq_len(M))
            rr3[[m]] <- outer(pp[, m], 1 - ff[, m])
        p3 <- Reduce("+", rr3) / M
        v3 <- apply(array(unlist(rr3), c(length(dt), nrow(cc), length(rr2))), c(1, 2), sd)^2 / M
        cov12 <- t(sapply(seq_along(dt), function (i) {
            cov(t(sapply(rr2, function (x) x[i, ])), pp[i, ])
        })) / M
        cov13 <- t(sapply(seq_along(dt), function (i) {
            cov(t(sapply(rr3, function (x) x[i, ])), pp[i, ])
        })) / M
        sn <- p2 / (1 - p1)
        sp <- p3 / p1
        se.sn <- sn
        se.sp <- sp
        for (k in 1:nrow(sn)) {
            dm <- function (x) {
                g <- c(x[1], x[2])
                V <- matrix(c(x[3], x[4], x[4], x[5]), 2, 2)
                suppressWarnings(t(g) %*% V %*% g)
            }
            ss <- cbind(1/(1-p1[k]), p2[k, ]/(1-p1[k])^2, v2[k, ], cov12[k, ], v1[k])
            se.sn[k, ] <- apply(ss, 1, dm)
            tt <- cbind(1/p1[k], -p3[k, ]/p1[k]^2, v3[k, ], cov13[k, ], v1[k])
            se.sp[k, ] <- apply(tt, 1, dm)
        }
        res[[i]] <- list("Sens" = sn, "Spec" = sp, "seSens" = se.sn, "seSpec" = se.sp)
    }
    aucs <- sapply(res, function (x) {
        sn <- t(x$Sen)
        fp <- t(1 - x$Spec)
        rr <- colSums(abs(diff(fp)) * sn[-1, ])
        names(rr) <- tail(sprintf("%.1f", dt), length(rr))
        rr
    })
    optThr <- lapply(res, function (x) {
        thr <- if (optThr == "sens*spec") t(x$Sen) * t(x$Spec) else t(x$Sen) + t(x$Spec) - 1
        ind <- apply(thr, 2, which.max)
        rr <- cc[ind, , drop = FALSE]
        rownames(rr) <- tail(sprintf("%.1f", dt), length(rr))
        rr
    })
    times <- split(data[[timeVar]], data[[idVar]][, drop = TRUE])
    for (i in seq_along(times))
        if (all(times[[i]] == 0))
            optThr[[i]] <- optThr[[i]][, 1, drop = FALSE]
    if (is.matrix(aucs)) {
        colnames(aucs) <- unique(data[[idVar]])
    } else {
        names(aucs) <- unique(data[[idVar]])
    }
    out <- c("MCresults" = res, list(AUCs = aucs, optThr = optThr, 
        times = times, dt = dt, M = M, diffType = diffType, 
        abs.diff = abs.diff, rel.diff = rel.diff, cc = cc, 
        min.cc = min.cc, max.cc = max.cc, success.rate = success.rate))
    class(out) <- "rocJM"
    out
}
