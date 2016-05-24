survfitJM.jointModel <- function (object, newdata, idVar = "id", simulate = TRUE, survTimes = NULL, 
            last.time = NULL, M = 200, CI.levels = c(0.025, 0.975), scale = 1.6, ...) {
    if (!inherits(object, "jointModel"))
        stop("Use only with 'jointModel' objects.\n")
    if (object$CompRisk)
        stop("survfitJM() is not currently implemented for ",
            "competing risks joint models.\n")
    if (!is.data.frame(newdata) || nrow(newdata) == 0)
        stop("'newdata' must be a data.frame with more than one rows.\n")
    if (is.null(newdata[[idVar]]))
        stop("'idVar' not in 'newdata.\n'")
    if (is.null(survTimes) || !is.numeric(survTimes))
        survTimes <- seq(min(exp(object$y$logT)), 
            max(exp(object$y$logT)) + 0.1, length.out = 35)
    method <- object$method
    timeVar <- object$timeVar
    interFact <- object$interFact
    parameterization <- object$parameterization
    derivForm <- object$derivForm
    indFixed <- derivForm$indFixed
    indRandom <- derivForm$indRandom
    LongFormat <- object$LongFormat
    TermsX <- object$termsYx
    TermsZ <- object$termsYz
    TermsX.deriv <- object$termsYx.deriv
    TermsZ.deriv <- object$termsYz.deriv
    mfX <- model.frame(TermsX, data = newdata)
    mfZ <- model.frame(TermsZ, data = newdata)
    formYx <- reformulate(attr(delete.response(TermsX), "term.labels"))
    formYz <- object$formYz
    na.ind <- as.vector(attr(mfX, "na.action"))
    na.ind <- if (is.null(na.ind)) {
        rep(TRUE, nrow(newdata))
    } else {
        !seq_len(nrow(newdata)) %in% na.ind
    }
    id <- as.numeric(unclass(newdata[[idVar]]))
    id <- id. <- match(id, unique(id))
    id <- id[na.ind]
    y <- model.response(mfX)
    X <- model.matrix(formYx, mfX)
    Z <- model.matrix(formYz, mfZ)[na.ind, , drop = FALSE]
    TermsT <- object$termsT
    data.id <- if (LongFormat) {
        nams.ind <- all.vars(delete.response(TermsT))
        ind <- !duplicated(newdata[nams.ind])
        newdata[ind, ]
    } else newdata[tapply(row.names(newdata), id, tail, n = 1L),]
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
    obs.times <- split(newdata[[timeVar]][na.ind], id)
    last.time <- if (is.null(last.time)) {
        tapply(newdata[[timeVar]], id., tail, n = 1)
    } else if (is.character(last.time) && length(last.time) == 1) {
        tapply(newdata[[last.time]], id., tail, n = 1)
    } else if (is.numeric(last.time) && length(last.time) == nrow(data.id)) {
        last.time
    } else {
        stop("\nnot appropriate value for 'last.time' argument.")
    }
    times.to.pred <- lapply(last.time, function (t) survTimes[survTimes > t])
    n <- object$n
    n.tp <- length(last.time)
    ncx <- ncol(X)
    ncz <- ncol(Z)
    ncww <- ncol(W)
    lag <- object$y$lag
    betas <- object$coefficients[['betas']]
    sigma <- object$coefficients[['sigma']]
    D <- object$coefficients[['D']]
    diag.D <- ncol(D) == 1 & nrow(D) > 1
    D <- if (diag.D) diag(c(D)) else D
    gammas <- object$coefficients[['gammas']]
    alpha <- object$coefficients[['alpha']]
    Dalpha <- object$coefficients[['Dalpha']]
    sigma.t <- object$coefficients[['sigma.t']]
    xi <- object$coefficients[['xi']]
    gammas.bs <- object$coefficients[['gammas.bs']]
    list.thetas <- list(betas = betas, log.sigma = log(sigma), gammas = gammas, alpha = alpha, 
            Dalpha = Dalpha, log.sigma.t = if (is.null(sigma.t)) NULL else log(sigma.t), 
            log.xi = if (is.null(xi)) NULL else log(xi), gammas.bs = gammas.bs, 
            D = if (diag.D) log(diag(D)) else chol.transf(D))
    if (method %in% c("weibull-PH-GH", "weibull-AFT-GH") && !is.null(object$scaleWB)) {
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
    if (!method %in% c("weibull-PH-GH", "weibull-AFT-GH", "piecewise-PH-GH", "spline-PH-GH")) {
        stop("\nsurvfitJM() is not yet available for this type of joint model.")
    }
    list.thetas <- list.thetas[!sapply(list.thetas, is.null)]
    thetas <- unlist(as.relistable(list.thetas))
    Var.thetas <- vcov(object)
    environment(log.posterior.b) <- environment(S.b) <- environment(ModelMats) <- environment()
    # construct model matrices to calculate the survival functions
    obs.times.surv <- split(data.id[[timeVar]], idT)
    survMats <- survMats.last <- vector("list", n.tp)
    for (i in seq_len(n.tp)) {
        survMats[[i]] <- lapply(times.to.pred[[i]], ModelMats, ii = i,
            obs.times = obs.times.surv, survTimes = survTimes)
        survMats.last[[i]] <- ModelMats(last.time[i], ii = i, 
            obs.times = obs.times.surv, survTimes = survTimes)
    }
    # calculate the Empirical Bayes estimates and their (scaled) variance
    modes.b <- matrix(0, n.tp, ncz)
    Vars.b <- vector("list", n.tp)
    for (i in seq_len(n.tp)) {
        betas.new <- betas
        sigma.new <- sigma
        D.new <- D
        gammas.new <- gammas
        alpha.new <- alpha
        Dalpha.new <- Dalpha
        sigma.t.new <- sigma.t
        xi.new <- xi
        gammas.bs.new <- gammas.bs
        ff <- function (b, y, tt, mm, i) -log.posterior.b(b, y, Mats = tt, method = mm, ii = i)
        opt <- try(optim(rep(0, ncz), ff, y = y, tt = survMats.last, mm = method, i = i, 
            method = "BFGS", hessian = TRUE), TRUE)
        if (inherits(opt, "try-error")) {
            gg <- function (b, y, tt, mm, i) cd(b, ff, y = y, tt = tt, mm = mm, i = i)
            opt <- optim(rep(0, ncz), ff, gg, y = y, tt = survMats.last, mm = method, 
                i = i, method = "BFGS", hessian = TRUE)
        } 
        modes.b[i, ] <- opt$par
        Vars.b[[i]] <- scale * solve(opt$hessian)        
    }
    if (!simulate) {
        res <- vector("list", n.tp)
        for (i in seq_len(n.tp)) {
            S.last <- S.b(last.time[i], modes.b[i, ], i, survMats.last[[i]])
            S.pred <- numeric(length(times.to.pred[[i]]))
            for (l in seq_along(S.pred))
                S.pred[l] <- S.b(times.to.pred[[i]][l], modes.b[i, ], i, survMats[[i]][[l]])
            res[[i]] <- cbind(times = times.to.pred[[i]], predSurv = S.pred / S.last)
            rownames(res[[i]]) <- seq_along(S.pred) 
        }
    } else {
        out <- vector("list", M)
        success.rate <- matrix(FALSE, M, n.tp)
        b.old <- b.new <- modes.b
        if (n.tp == 1)
            dim(b.old) <- dim(b.new) <- c(1, ncz)    
        for (m in 1:M) {
            # Step 1: simulate new parameter values
            thetas.new <- mvrnorm(1, thetas, Var.thetas)
            thetas.new <- relist(thetas.new, skeleton = list.thetas)
            betas.new <- thetas.new$betas
            sigma.new <- exp(thetas.new$log.sigma)
            gammas.new <- thetas.new$gammas
            alpha.new <- thetas.new$alpha
            D.new <- thetas.new$D
            D.new <- if (diag.D) exp(D.new) else chol.transf(D.new)
            if (method == "weibull-PH-GH" || method == "weibull-AFT-GH") {
                sigma.t.new <- if (is.null(object$scaleWB)) exp(thetas.new$log.sigma.t) else object$scaleWB
            } else if (method == "piecewise-PH-GH") {
                xi.new <- exp(thetas.new$log.xi)
            } else if (method == "spline-PH-GH") {
                gammas.bs.new <- thetas.new$gammas.bs
            }
            SS <- vector("list", n.tp)
            for (i in seq_len(n.tp)) {
                # Step 2: simulate new random effects values
                proposed.b <- rmvt(1, modes.b[i, ], Vars.b[[i]], 4)
                dmvt.old <- dmvt(b.old[i, ], modes.b[i, ], Vars.b[[i]], 4, TRUE)
                dmvt.proposed <- dmvt(proposed.b, modes.b[i, ], Vars.b[[i]], 4, TRUE)
                a <- min(exp(log.posterior.b(proposed.b, y, survMats.last, method, ii = i) + dmvt.old - 
                        log.posterior.b(b.old[i, ], y, survMats.last, method, ii = i) - dmvt.proposed), 1)
                ind <- runif(1) <= a
                success.rate[m, i] <- ind
                if (!is.na(ind) && ind)
                    b.new[i, ] <- proposed.b
                # Step 3: compute Pr(T > t_k | T > t_{k - 1}; theta.new, b.new)
                S.last <- S.b(last.time[i], b.new[i, ], i, survMats.last[[i]])
                S.pred <- numeric(length(times.to.pred[[i]]))
                for (l in seq_along(S.pred))
                    S.pred[l] <- S.b(times.to.pred[[i]][l], b.new[i, ], i, survMats[[i]][[l]])
                SS[[i]] <- S.pred / S.last
            }
            b.old <- b.new
            out[[m]] <- SS
        }
        res <- vector("list", n.tp)
        for (i in seq_len(n.tp)) {
            rr <- sapply(out, "[[", i)
            if (!is.matrix(rr))
                rr <- rbind(rr)
            res[[i]] <- cbind(
                times = times.to.pred[[i]],
                "Mean" = rowMeans(rr, na.rm = TRUE),
                "Median" = apply(rr, 1, median, na.rm = TRUE),
                "Lower" = apply(rr, 1, quantile, probs = CI.levels[1], na.rm = TRUE),
                "Upper" = apply(rr, 1, quantile, probs = CI.levels[2], na.rm = TRUE)
            )
            rownames(res[[i]]) <- as.character(seq_len(NROW(res[[i]])))
        }
    }
    y <- split(y, id)
    newdata. <- do.call(rbind, mapply(function (d, t) {
        d. <- rbind(d, d[nrow(d), ])
        d.[[timeVar]][nrow(d.)] <- t
        d.
    }, split(newdata, id.), last.time, SIMPLIFY = FALSE))
    id. <- as.numeric(unclass(newdata.[[idVar]]))
    id. <- match(id., unique(id.))
    mfX. <- model.frame(delete.response(TermsX), data = newdata.)
    mfZ. <- model.frame(TermsZ, data = newdata.)
    X. <- model.matrix(formYx, mfX.)
    Z. <- model.matrix(formYz, mfZ.)
    fitted.y <- split(c(X. %*% betas) + rowSums(Z. * modes.b[id., , drop = FALSE]), id.)
    names(res) <- names(y) <- names(last.time) <- names(obs.times) <- unique(unclass(newdata[[idVar]]))
    res <- list(summaries = res, survTimes = survTimes, last.time = last.time, 
        obs.times = obs.times, y = y, 
        fitted.times = split(newdata.[[timeVar]], factor(newdata.[[idVar]])), 
        fitted.y = fitted.y, ry = range(object$y$y, na.rm = TRUE))
    if (simulate) {
        res$full.results <- out
        res$success.rate <- success.rate
    }
    rm(list = ".Random.seed", envir = globalenv())
    class(res) <- "survfitJM"
    res
}
