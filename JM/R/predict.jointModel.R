predict.jointModel <-
function (object, newdata, type = c("Marginal", "Subject"),
    interval = c("none", "confidence", "prediction"), level = 0.95, idVar = "id", 
    FtTimes = NULL, M = 300, returnData = FALSE, scale = 1.6, ...) {
    if (!inherits(object, "jointModel"))
        stop("Use only with 'jointModel' objects.\n")
    if (!is.data.frame(newdata) || nrow(newdata) == 0)
        stop("'newdata' must be a data.frame with more than one rows.\n")
    type <- match.arg(type)
    interval <- match.arg(interval)
    if (type == "Marginal") {
        TermsX <- delete.response(object$termsYx)
        mf <- model.frame(TermsX, data = newdata)
        form <- reformulate(attr(TermsX, "term.labels"))
        X <- model.matrix(form, data = mf)
        out <- c(X %*% object$coefficients$betas)
        names(out) <- row.names(newdata)
        if (interval == "prediction") {
            warning("\nfor type = 'Marginal' only confidence intervals are calculated.")
            interval <- "confidence"
        }
        if (interval == "confidence") {
            V <- vcov(object)
            ind <- head(grep("Y.", colnames(V), fixed = TRUE), -1)
            se.fit <- sqrt(diag(X %*% tcrossprod(V[ind, ind], X)))
            alpha <- 1 - level
            low <- out + qnorm(alpha/2) * se.fit
            up <- out + qnorm(1-alpha/2) * se.fit
            names(se.fit) <- names(low) <- names(up) <- row.names(newdata)
            out <- list(pred = out, se.fit = se.fit, low = low, upp = up)
        }
        if (returnData) {
            out <- if (is.list(out)) 
                cbind(newdata, do.call(cbind, out))
            else
                cbind(newdata, pred = out)
        }
    } else {
        if (object$CompRisk)
            stop("predict() with type = 'Subject' is not currently ",
                "\n  implemented for competing risks joint models.\n")
        if (is.null(newdata[[idVar]]))
            stop("'idVar' not in 'newdata.\n'")
        method <- object$method
        timeVar <- object$timeVar
        interFact <- object$interFact
        parameterization <- object$parameterization
        derivForm <- object$derivForm
        indFixed <- derivForm$indFixed
        indRandom <- derivForm$indRandom
        #id <- as.numeric(unclass(newdata[[idVar]]))
        #id <- match(id, unique(id))
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
        #data.id <- if (LongFormat) newdata else newdata[!duplicated(id), ]
        data.id <- if (LongFormat) {
            nams.ind <- all.vars(delete.response(TermsT))
            ind <- !duplicated(newdata[nams.ind])
            newdata[ind, ]
        } else newdata[!duplicated(id), ]
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
        obs.times <- split(newdata[[timeVar]], id.)
        last.time <- tapply(newdata[[timeVar]], id., tail, n = 1)
        times.to.pred <- if (is.null(FtTimes)) {
            lapply(last.time, 
                function (t) seq(t, max(object$times) + 
                    0.1 * mad(object$times), length = 25))
        } else {
            if (!is.list(FtTimes) || length(FtTimes) != length(last.time))
                rep(list(FtTimes), length(last.time))
            else
                FtTimes
        }
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
        list.thetas <- list(betas = betas, log.sigma = log(sigma), 
            gammas = gammas, alpha = alpha, 
            Dalpha = Dalpha, 
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
            stop("\nsurvfitJM() is not yet available for this type of joint model.")
        }
        list.thetas <- list.thetas[!sapply(list.thetas, is.null)]
        thetas <- unlist(as.relistable(list.thetas))
        Var.thetas <- vcov(object)
        environment(log.posterior.b) <- environment(ModelMats) <- environment()
        # construct model matrices to calculate the survival functions
        obs.times.surv <- split(data.id[[timeVar]], idT)
        survMats.last <- vector("list", n.tp)
        for (i in seq_len(n.tp)) {
            survMats.last[[i]] <- ModelMats(last.time[i], ii = i,
                obs.times = obs.times.surv, 
                survTimes = unlist(times.to.pred, use.names = FALSE))
        }
        data.id2 <- newdata[!duplicated(id), ]
        data.id2 <- data.id2[rep(1:nrow(data.id2), 
            sapply(times.to.pred, length)), ]
        data.id2[[timeVar]] <- unlist(times.to.pred)
        mfXpred <- model.frame(TermsX, data = data.id2)
        mfZpred <- model.frame(TermsZ, data = data.id2)
        Xpred <- model.matrix(formYx, mfXpred)
        Zpred <- model.matrix(formYz, mfZpred)
        id2 <- as.numeric(unclass(data.id2[[idVar]]))
        id2 <- match(id2, unique(id2))
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
            ff <- function (b, y, tt, mm, i) 
                -log.posterior.b(b, y, Mats = tt, method = mm, ii = i)
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
        res <- vector("list", M)
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
                sigma.t.new <- if (is.null(object$scaleWB)) 
                    exp(thetas.new$log.sigma.t) else object$scaleWB
            } else if (method == "piecewise-PH-GH") {
                xi.new <- exp(thetas.new$log.xi)
            } else if (method == "spline-PH-GH") {
                gammas.bs.new <- thetas.new$gammas.bs
            }
            y.new <- vector("list", n.tp)
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
                # Step 3: compute future Ys
                Xpred.i <- Xpred[id2 == i, , drop = FALSE]
                Zpred.i <- Zpred[id2 == i, , drop = FALSE]
                mu.i <- as.vector(c(Xpred.i %*% betas.new) + 
                    rowSums(Zpred.i * rep(b.new[i, ], each = nrow(Zpred.i))))
                y.new[[i]] <- if (interval == "confidence") mu.i else 
                    if (interval == "prediction") rnorm(length(mu.i), mu.i, sigma.new)
            }
            b.old <- b.new
            res[[m]] <- y.new
        }
        oo <- vector("list", n.tp)
        for (i in seq_len(n.tp)) {
            oo[[i]] <- do.call(rbind, sapply(res, "[", i))
        }
        out <- as.vector(c(Xpred %*% betas) + 
            rowSums(Zpred * modes.b[id2, , drop = FALSE]))
        if (interval %in% c("confidence", "prediction")) {
            alpha <- 1 - level
            se.fit <- lapply(oo, function (m) {
                if (is.matrix(m)) 
                    apply(m, 2, sd)
                else 
                    sd(m)
            })
            f1 <- function (mat) apply(mat, 2, quantile, probs = alpha/2)
            f2 <- function (mat) apply(mat, 2, quantile, probs = 1 - alpha/2)
            low <- lapply(oo, f1) 
            up <- lapply(oo, f2)
            out <- list(pred = out, se.fit = unlist(se.fit), 
                low = unlist(low), upp = unlist(up))
        }
        if (returnData) {
            newdata$pred <- c(X %*% betas) + rowSums(Z * modes.b[id, ])
            out <- if (is.list(out)) {
                newdata$upp <- newdata$low <- newdata$se.fit <- NA
                rbind(newdata, cbind(data.id2, do.call(cbind, out)))
            } else {
                rbind(newdata, cbind(data.id2, pred = out))
            }
        } else
            attr(out, "time.to.pred") <- times.to.pred
    }
    rm(list = ".Random.seed", envir = globalenv())
    out
}
