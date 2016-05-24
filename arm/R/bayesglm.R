bayesglm <- function (formula, family = gaussian, data, weights, subset,
    na.action, start = NULL, etastart, mustart, offset, control = list(...),
    model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, contrasts = NULL,
    drop.unused.levels = TRUE,
    prior.mean = 0, prior.scale = NULL, prior.df = 1,
    prior.mean.for.intercept = 0,
    prior.scale.for.intercept = NULL,
    prior.df.for.intercept = 1,
    min.prior.scale = 1e-12,
    scaled = TRUE,
    keep.order = TRUE, drop.baseline = TRUE,
    maxit = 100, print.unnormalized.log.posterior = FALSE,
    Warning = TRUE, ...)
{
    call <- match.call()
    if (is.character(family)) {
        family <- get(family, mode = "function", envir = parent.frame())
    }
    if (is.function(family)) {
        family <- family()
    }
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
    if (missing(data)) {
        data <- environment(formula)
    }
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action",
        "etastart", "mustart", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- drop.unused.levels
    mf$na.action <- NULL
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    if (identical(method, "model.frame")){
        return(mf)
    }
    if (!is.character(method) && !is.function(method)){
        stop("invalid 'method' argument")
    }
    if (identical(method, "glm.fit")){
        control <- do.call("glm.control", control)
    }
    control$maxit <- maxit

    mt <- attr(mf, "terms")
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm)) {
            names(Y) <- nm
        }
    }
    X <- if (!is.empty.model(mt)) {
        model.matrixBayes(object = mt, data = data, contrasts.arg = contrasts,
            keep.order = keep.order, drop.baseline = drop.baseline)
    }else {
        matrix(, NROW(Y), 0L)
    }
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights)) {
        stop("'weights' must be a numeric vector")
    }
    if (!is.null(weights) && any(weights < 0)) {
        stop("negative weights not allowed")
    }
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(Y))
            stop(gettextf("number of offsets is %d should equal %d (number of
            observations)",
                length(offset), NROW(Y)), domain = NA)
    }
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")
    fit <- bayesglm.fit(x = X, y = Y, weights = weights, start = start,
        etastart = etastart, mustart = mustart, offset = offset,
        family = family, control = control,
        intercept = attr(mt, "intercept") > 0L,
        prior.mean = prior.mean,
        prior.scale = prior.scale,
        prior.df = prior.df,
        prior.mean.for.intercept = prior.mean.for.intercept,
        prior.scale.for.intercept = prior.scale.for.intercept,
        prior.df.for.intercept = prior.df.for.intercept,
        min.prior.scale = min.prior.scale,
        print.unnormalized.log.posterior = print.unnormalized.log.posterior,
        scaled = scaled, Warning = Warning)
    if (length(offset) && attr(mt, "intercept") > 0L) {
        fit2 <- bayesglm.fit(x = X[, "(Intercept)",
            drop = FALSE], y = Y, weights = weights, offset = offset,
            family = family, control = control, intercept = TRUE,
            prior.mean = prior.mean, prior.scale = prior.scale,
            prior.df = prior.df, prior.mean.for.intercept =
            prior.mean.for.intercept,
            prior.scale.for.intercept = prior.scale.for.intercept,
            prior.df.for.intercept = prior.df.for.intercept,
            min.prior.scale = min.prior.scale,
            print.unnormalized.log.posterior = print.unnormalized.log.posterior,
            scaled = scaled,
            Warning = Warning)
        if (!fit2$converged){
            warning("fitting to calculate the null deviance did not converge --
            increase 'maxit'?")
        }
        fit$null.deviance <- fit2$deviance
    }
    if (model) {
        fit$model <- mf
    }
    fit$na.action <- attr(mf, "na.action")
    if (x) {
        fit$x <- X
    }
    if (!y) {
        fit$y <- NULL
    }
    fit <- c(fit, list(call = call, formula = formula, terms = mt,
        data = data, offset = offset, control = control, method = method,
        contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt,
            mf)), keep.order = keep.order, drop.baseline = drop.baseline)
    class(fit) <- c("bayesglm", "glm", "lm")
    return(fit)
}



bayesglm.fit <- function (x, y, weights = rep(1, nobs), start = NULL, etastart = NULL,
    mustart = NULL, offset = rep(0, nobs), family = gaussian(),
    control = list(), intercept = TRUE,
    prior.mean = 0, prior.scale = NULL, prior.df = 1,
    prior.mean.for.intercept = 0,
    prior.scale.for.intercept = NULL,
    prior.df.for.intercept = 1,
    min.prior.scale = 1e-12, scaled = TRUE,
    print.unnormalized.log.posterior = FALSE,
    Warning = TRUE)
{

    control <- do.call("glm.control", control)

    x <- as.matrix(x)
    xnames <- dimnames(x)[[2L]]
    ynames <- if (is.matrix(y)){
        rownames(y)
    }else{
        names(y)
    }
    conv <- FALSE
    nobs <- NROW(y)
    nvars <- NCOL(x)
#===============================
#  initialize priors
#===============================
    if(is.null(prior.scale)){
        prior.scale <- 2.5
        if(family$link == "probit"){
            prior.scale <- prior.scale*1.6
        }
    }

    if(is.null(prior.scale.for.intercept)){
        prior.scale.for.intercept <- 10
        if(family$link == "probit"){
            prior.scale.for.intercept <- prior.scale.for.intercept*1.6
        }
    }

    if(intercept){
        nvars <- nvars - 1
    }

    if(length(prior.mean)==1L){
        prior.mean <- rep(prior.mean, nvars)
    }else if(length(prior.mean)!=nvars){
        stop("invalid length for prior.mean")
    }

    if(length(prior.scale)==1L){
        prior.scale <- rep(prior.scale, nvars)
    }else if(length(prior.scale)!=nvars){
        stop("invalid length for prior.scale")
    }

    if(length(prior.df)==1L){
        prior.df <- rep(prior.df, nvars)
    }else if(length(prior.df)!=nvars){
        stop("invalid length for prior.df")
    }

    if(intercept){
        prior.mean <- c(prior.mean.for.intercept, prior.mean)
        prior.scale <- c(prior.scale.for.intercept, prior.scale)
        prior.df <- c(prior.df.for.intercept, prior.df)
    }

    if(scaled){
        if(family$family=="gaussian"){
            prior.scale <- prior.scale*2*sd(y)
        }
        prior.scale.0 <- prior.scale
        if(nvars==0) nvars = 1
        for(j in 1:nvars){
            x.obs <- x[,j]
            x.obs <- x.obs[!is.na(x.obs)]
            num.categories <- length(unique(x.obs))
            x.scale <- 1
            if(num.categories==2L){
                x.scale <- max(x.obs) - min(x.obs)
            }else if(num.categories>2){
                x.scale <- 2*sd(x.obs)
            }
            prior.scale[j] <- prior.scale[j]/x.scale
            if(prior.scale[j] < min.prior.scale){
                prior.scale[j] <- min.prior.scale
                warning("prior scale for varible ", j,
                    " set to min.prior.scale = ", min.prior.scale, "\n")
            }
        }
    }

#===================
    nvars <- NCOL(x)
    EMPTY <- nvars == 0
    if (is.null(weights))
        weights <- rep.int(1, nobs)
    if (is.null(offset))
        offset <- rep.int(0, nobs)
    variance <- family$variance
    linkinv <- family$linkinv
    if (!is.function(variance) || !is.function(linkinv))
        stop("'family' argument seems not to be a valid family object",
            call. = FALSE)
    dev.resids <- family$dev.resids
    aic <- family$aic
    mu.eta <- family$mu.eta
    unless.null <- function(x, if.null){
        if (is.null(x))
            if.null
        else x
    }
    valideta <- unless.null(family$valideta, function(eta) TRUE)
    validmu <- unless.null(family$validmu, function(mu) TRUE)
    if (is.null(mustart)) {
        eval(family$initialize)
    }else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }
    if (EMPTY) {
        eta <- rep.int(0, nobs) + offset
        if (!valideta(eta))
            stop("invalid linear predictor values in empty model",
                call. = FALSE)
        mu <- linkinv(eta)
        if (!validmu(mu))
            stop("invalid fitted means in empty model", call. = FALSE)
        dev <- sum(dev.resids(y, mu, weights))
        w <- ((weights * mu.eta(eta)^2)/variance(mu))^0.5
        residuals <- (y - mu)/mu.eta(eta)
        good <- rep_len(TRUE, length(residuals))
        boundary <- conv <- TRUE
        coef <- numeric()
        iter <- 0L
    } else {
        coefold <- NULL
        eta <- if (!is.null(etastart)){
            etastart
        }else if (!is.null(start)){
            if (length(start) != nvars){
                if(start==0&length(start)==1){
                    start <- rep(0, nvars)
                    offset + as.vector(ifelse((NCOL(x) == 1L), x*start, x %*% start))
                }else{
                    stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s",
                    nvars, paste(deparse(xnames), collapse = ", ")),
                    domain = NA)
                }
            } else {
                coefold <- start
                offset + as.vector(if (NCOL(x) == 1L)
                  x * start
                else x %*% start)
            }
        }else{
            family$linkfun(mustart)
        }
        mu <- linkinv(eta)
        if (!(validmu(mu) && valideta(eta)))
            stop("cannot find valid starting values: please specify some",
                call. = FALSE)
        devold <- sum(dev.resids(y, mu, weights))
        boundary <- conv <- FALSE

#======================================
#   initialize prior.sd
#======================================
    prior.sd <- prior.scale
#=====================================
    dispersion <- ifelse((family$family %in% c("poisson", "binomial")),  1, var(y)/10000)
    dispersionold <- dispersion
    for (iter in 1L:control$maxit) {
        good <- weights > 0
        varmu <- variance(mu)[good]
        if (anyNA(varmu))
            stop("NAs in V(mu)")
        if (any(varmu == 0))
            stop("0s in V(mu)")
        mu.eta.val <- mu.eta(eta)
        if (any(is.na(mu.eta.val[good])))
            stop("NAs in d(mu)/d(eta)")
        good <- (weights > 0) & (mu.eta.val != 0)
        if (all(!good)) {
            conv <- FALSE
            warning(gettextf("no observations informative at iteration %d",
                iter), domain = NA)
                break
        }
        z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
        w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
        ngoodobs <- as.integer(nobs - sum(!good))
#======================
#  data augmentation
#=========================
        # coefs.hat <- rep(0, NCOL(x))   # why do we need coefs.hat here? SU 2015.3.30
        x.star <- rbind(x, diag(NCOL(x)))
        if(intercept&scaled){
            x.star[nobs+1,] <- colMeans(x)
        }
        z.star <- c (z, prior.mean)
        w.star <- c (w, sqrt(dispersion)/prior.sd)
#=================================================
        good.star <- c (good, rep(TRUE,NCOL(x)))
        ngoodobs.star <- ngoodobs + NCOL(x)
            #fit <- .Call(C_Cdqrls, x.star[good, , drop = FALSE] *
            #    w.star, z.star * w.star, min(1e-07, control$epsilon/1000),
            #    check = FALSE)
        fit <- lm.fit(x = x.star[good.star,,drop=FALSE]*w.star, y = z.star*w.star)
        if (any(!is.finite(fit$coefficients))) {
            conv <- FALSE
            warning(gettextf("non-finite coefficients at iteration %d",
                  iter), domain = NA)
            break
        }
        start[fit$qr$pivot] <- coefs.hat <- fit$coefficients
        fit$qr$qr <- as.matrix (fit$qr$qr)
        V.coefs <- chol2inv(fit$qr$qr[1:NCOL(x.star), 1:NCOL(x.star), drop = FALSE])
        if (family$family == "gaussian" & scaled){
            prior.scale <- prior.scale.0
        }
        prior.sd <- ifelse(prior.df == Inf, prior.scale,
            sqrt(((coefs.hat - prior.mean)^2 + diag(V.coefs)*dispersion +
            prior.df * prior.scale^2)/(1 + prior.df)))
        start[fit$qr$pivot] <- fit$coefficients
        eta <- drop(x %*% start)
        mu <- linkinv(eta <- eta + offset)
        dev <- sum(dev.resids(y, mu, weights))
        if (!(family$family %in% c("poisson", "binomial"))) {
            mse.resid <- mean((w * (z - x %*% coefs.hat))^2)
            mse.uncertainty <- mean(rowSums(( x %*% V.coefs ) * x)) * dispersion # faster
            dispersion <- mse.resid + mse.uncertainty
        }
        if (control$trace)
            cat("Deviance = ", dev, " Iterations - ", iter,
                "\n", sep = "")
        boundary <- FALSE
        if (!is.finite(dev)) {
            if (is.null(coefold))
                stop("no valid set of coefficients has been found: please supply starting values",
                    call. = FALSE)
                warning("step size truncated due to divergence",
                  call. = FALSE)
                ii <- 1
                while (!is.finite(dev)) {
                  if (ii > control$maxit)
                    stop("inner loop 1; cannot correct step size",
                      call. = FALSE)
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- drop(x %*% start)
                  mu <- linkinv(eta <- eta + offset)
                  dev <- sum(dev.resids(y, mu, weights))
                }
                boundary <- TRUE
                if (control$trace)
                  cat("Step halved: new deviance = ", dev, "\n",
                    sep = "")
        }
        if (!(valideta(eta) && validmu(mu))) {
                if (is.null(coefold))
                  stop("no valid set of coefficients has been found: please supply starting values",
                    call. = FALSE)
                warning("step size truncated: out of bounds",
                  call. = FALSE)
                ii <- 1
                while (!(valideta(eta) && validmu(mu))) {
                  if (ii > control$maxit)
                    stop("inner loop 2; cannot correct step size",
                      call. = FALSE)
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- drop(x %*% start)
                  mu <- linkinv(eta <- eta + offset)
                }
                boundary <- TRUE
                dev <- sum(dev.resids(y, mu, weights))
                if (control$trace)
                  cat("Step halved: new deviance = ", dev, "\n",
                    sep = "")
        }
#===============================
# print unnormalized log posterior
#================================
        if (family$family == "binomial" && print.unnormalized.log.posterior) {
            logprior <- sum(dt(coefs.hat, prior.df, prior.mean, log = TRUE))
            xb <- invlogit( x %*% coefs.hat )
            loglikelihood <- sum( log( c( xb[ y == 1 ], 1 - xb[ y == 0 ] ) ) )
            cat( "log prior: ", logprior, ", log likelihood: ", loglikelihood, ",
            unnormalized log posterior: ", loglikelihood +logprior, "\n" ,sep="")
        }
#================================

        if (iter > 1 & abs(dev - devold)/(0.1 + abs(dev)) <
                control$epsilon & abs(dispersion - dispersionold)/(0.1 +
                abs(dispersion)) < control$epsilon) {
                conv <- TRUE
                coef <- start
                break
        }else {
            devold <- dev
            dispersionold <- dispersion
            coef <- coefold <- start
        }
    }
        if (!conv){
            warning("algorithm did not converge", call. = FALSE)
        }
        if (boundary){
            warning("algorithm stopped at boundary value",
                call. = FALSE)
        }
        eps <- 10 * .Machine$double.eps
        if (family$family == "binomial") {
            if (any(mu > 1 - eps) || any(mu < eps)){
                warning("fitted probabilities numerically 0 or 1 occurred",
                  call. = FALSE)
            }
        }
        if (family$family == "poisson") {
            if (any(mu < eps)){
                warning("fitted rates numerically 0 occurred",
                  call. = FALSE)
            }
        }
        if (fit$rank < nvars){
            coef[fit$qr$pivot][seq.int(fit$rank + 1, nvars)] <- NA
        }
        xxnames <- xnames[fit$qr$pivot]
        residuals <- rep.int(NA, nobs)
        residuals[good] <- z - (eta - offset)[good]
        fit$qr$qr <- as.matrix(fit$qr$qr)
        nr <- min(sum(good), nvars)
        if (nr < nvars) {
            Rmat <- diag(nvars)
            Rmat[1L:nr, 1L:nvars] <- fit$qr$qr[1L:nr, 1L:nvars]
        } else Rmat <- fit$qr$qr[1L:nvars, 1L:nvars]
        Rmat <- as.matrix(Rmat)
        Rmat[row(Rmat) > col(Rmat)] <- 0
        names(coef) <- xnames
        colnames(fit$qr$qr) <- xxnames
        dimnames(Rmat) <- list(xxnames, xxnames)
    }
    names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    wt <- rep.int(0, nobs)
    wt[good] <- w^2
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
    wtdmu <- if (intercept){
        sum(weights * y)/sum(weights)
    } else{
        linkinv(offset)
    }
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(intercept)
    rank <- if (EMPTY) {
        0
    } else{
        fit$rank
    }
    resdf <- n.ok - rank
    aic.model <- aic(y, n.ok, mu, weights, dev) + 2 * rank
    list(coefficients = coef,
        residuals = residuals,
        fitted.values = mu,
        effects = if (!EMPTY) fit$effects,
        R = if (!EMPTY) Rmat,
        rank = rank,
        qr = if (!EMPTY) structure(getQr(fit)[c("qr", "rank", "qraux", "pivot", "tol")], class = "qr"),
        family = family,
        linear.predictors = eta,
        deviance = dev,
        aic = aic.model,
        null.deviance = nulldev,
        iter = iter,
        weights = wt,
        prior.weights = weights,
        df.residual = resdf,
        df.null = nulldf,
        y = y,
        converged = conv,
        boundary = boundary,
        prior.mean = prior.mean,
        prior.scale = prior.scale,
        prior.df = prior.df,
        prior.sd = prior.sd,
        dispersion = dispersion)
}

setMethod("print", signature(x = "bayesglm"),
    function(x, digits=2) display(object=x, digits=digits))

setMethod("show", signature(object = "bayesglm"),
    function(object) display(object, digits=2))

predict.bayesglm <- function (object, newdata = NULL, type = c("link", "response",
    "terms"), se.fit = FALSE, dispersion = NULL, terms = NULL,
    na.action = na.pass, ...)
{
    type <- match.arg(type)
    na.act <- object$na.action
    object$na.action <- NULL
    if (!se.fit) {
        if (missing(newdata)) {
            pred <- switch(type, link = object$linear.predictors,
                response = object$fitted.values, terms = predictLM(object,
                  se.fit = se.fit, scale = 1, type = "terms",
                  terms = terms))
            if (!is.null(na.act))
                pred <- napredict(na.act, pred)
        }
        else {
            pred <- predictLM(object, newdata, se.fit, scale = 1,
                type = ifelse(type == "link", "response", type),
                terms = terms, na.action = na.action)
            switch(type, response = {
                pred <- family(object)$linkinv(pred)
            }, link = , terms = )
        }
    }
    else {
        if (inherits(object, "survreg"))
            dispersion <- 1
        if (is.null(dispersion) || dispersion == 0)
            dispersion <- summary(object, dispersion = dispersion)$dispersion
        residual.scale <- as.vector(sqrt(dispersion))
        pred <- predictLM(object, newdata, se.fit, scale = residual.scale,
            type = ifelse(type == "link", "response", type),
            terms = terms, na.action = na.action)
        fit <- pred$fit
        se.fit <- pred$se.fit
        switch(type, response = {
            se.fit <- se.fit * abs(family(object)$mu.eta(fit))
            fit <- family(object)$linkinv(fit)
        }, link = , terms = )
        if (missing(newdata) && !is.null(na.act)) {
            fit <- napredict(na.act, fit)
            se.fit <- napredict(na.act, se.fit)
        }
        pred <- list(fit = fit, se.fit = se.fit, residual.scale = residual.scale)
    }
    pred
}

predictLM <- function (object, newdata, se.fit = FALSE, scale = NULL, df = Inf,
    interval = c("none", "confidence", "prediction"), level = 0.95,
    type = c("response", "terms"), terms = NULL, na.action = na.pass,
    pred.var = res.var/weights, weights = 1, ...)
{
    tt <- terms(object)
    keep.order <- object$keep.order
    drop.baseline <- object$drop.baseline
    if (!inherits(object, "lm"))
        warning("calling predict.lm(<fake-lm-object>) ...")
    if (missing(newdata) || is.null(newdata)) {
        mm <- X <- model.matrix(object)
        mmDone <- TRUE
        offset <- object$offset
    }
    else {
        Terms <- delete.response(tt)
        m <- model.frame(Terms, newdata, na.action = na.action,
            xlev = object$xlevels)
        if (!is.null(cl <- attr(Terms, "dataClasses")))
            .checkMFClasses(cl, m)
        X <- model.matrixBayes(Terms, m, contrasts.arg = object$contrasts, keep.order = keep.order, drop.baseline = drop.baseline)
        offset <- rep(0, nrow(X))
        if (!is.null(off.num <- attr(tt, "offset")))
            for (i in off.num) offset <- offset + eval(attr(tt,
                "variables")[[i + 1]], newdata)
        if (!is.null(object$call$offset))
            offset <- offset + eval(object$call$offset, newdata)
        mmDone <- FALSE
    }
    n <- length(object$residuals)
    p <- object$rank
    p1 <- seq_len(p)
    piv <- if (p)
        getQr(object)$pivot[p1]
    if (p < ncol(X) && !(missing(newdata) || is.null(newdata)))
        warning("prediction from a rank-deficient fit may be misleading")
    beta <- object$coefficients
    predictor <- drop(X[, piv, drop = FALSE] %*% beta[piv])
    if (!is.null(offset))
        predictor <- predictor + offset
    interval <- match.arg(interval)
    if (interval == "prediction") {
        if (missing(newdata))
            warning("Predictions on current data refer to _future_ responses\n")
        if (missing(newdata) && missing(weights)) {
            w <- .weights.default(object)
            if (!is.null(w)) {
                weights <- w
                warning("Assuming prediction variance inversely proportional to weights used for fitting\n")
            }
        }
        if (!missing(newdata) && missing(weights) && !is.null(object$weights) &&
            missing(pred.var))
            warning("Assuming constant prediction variance even though model fit is weighted\n")
        if (inherits(weights, "formula")) {
            if (length(weights) != 2L)
                stop("'weights' as formula should be one-sided")
            d <- if (missing(newdata) || is.null(newdata))
                model.frame(object)
            else newdata
            weights <- eval(weights[[2L]], d, environment(weights))
        }
    }
    type <- match.arg(type)
    if (se.fit || interval != "none") {
        res.var <- if (is.null(scale)) {
            r <- object$residuals
            w <- object$weights
            rss <- sum(if (is.null(w)) r^2 else r^2 * w)
            df <- object$df.residual
            rss/df
        }
        else scale^2
        if (type != "terms") {
            if (p > 0) {
                XRinv <- if (missing(newdata) && is.null(w))
                  qr.Q(getQr(object))[, p1, drop = FALSE]
                else X[, piv] %*% qr.solve(qr.R(getQr(object))[p1, p1])
                ip <- drop(XRinv^2 %*% rep(res.var, p))
            }
            else ip <- rep(0, n)
        }
    }
    if (type == "terms") {
        if (!mmDone) {
            mm <- model.matrixBayes(object, keep.order = keep.order, drop.baseline = drop.baseline)
            mmDone <- TRUE
        }
        aa <- attr(mm, "assign")
        ll <- attr(tt, "term.labels")
        hasintercept <- attr(tt, "intercept") > 0L
        if (hasintercept)
            ll <- c("(Intercept)", ll)
        aaa <- factor(aa, labels = ll)
        asgn <- split(order(aa), aaa)
        if (hasintercept) {
            asgn$"(Intercept)" <- NULL
            if (!mmDone) {
                mm <- model.matrixBayes(object, keep.order = keep.order, drop.baseline = drop.baseline)
                mmDone <- TRUE
            }
            avx <- colMeans(mm)
            termsconst <- sum(avx[piv] * beta[piv])
        }
        nterms <- length(asgn)
        if (nterms > 0) {
            predictor <- matrix(ncol = nterms, nrow = NROW(X))
            dimnames(predictor) <- list(rownames(X), names(asgn))
            if (se.fit || interval != "none") {
                ip <- matrix(ncol = nterms, nrow = NROW(X))
                dimnames(ip) <- list(rownames(X), names(asgn))
                Rinv <- qr.solve(qr.R(getQr(object))[p1, p1])
            }
            if (hasintercept)
                X <- sweep(X, 2L, avx, check.margin = FALSE)
            unpiv <- rep.int(0L, NCOL(X))
            unpiv[piv] <- p1
            for (i in seq.int(1L, nterms, length.out = nterms)) {
                iipiv <- asgn[[i]]
                ii <- unpiv[iipiv]
                iipiv[ii == 0L] <- 0L
                predictor[, i] <- if (any(iipiv > 0L))
                  X[, iipiv, drop = FALSE] %*% beta[iipiv]
                else 0
                if (se.fit || interval != "none")
                  ip[, i] <- if (any(iipiv > 0L))
                    as.matrix(X[, iipiv, drop = FALSE] %*% Rinv[ii,
                      , drop = FALSE])^2 %*% rep.int(res.var,
                      p)
                  else 0
            }
            if (!is.null(terms)) {
                predictor <- predictor[, terms, drop = FALSE]
                if (se.fit)
                  ip <- ip[, terms, drop = FALSE]
            }
        }
        else {
            predictor <- ip <- matrix(0, n, 0L)
        }
        attr(predictor, "constant") <- if (hasintercept)
            termsconst
        else 0
    }
    if (interval != "none") {
        tfrac <- qt((1 - level)/2, df)
        hwid <- tfrac * switch(interval, confidence = sqrt(ip),
            prediction = sqrt(ip + pred.var))
        if (type != "terms") {
            predictor <- cbind(predictor, predictor + hwid %o%
                c(1, -1))
            colnames(predictor) <- c("fit", "lwr", "upr")
        }
        else {
            if (!is.null(terms))
                hwid <- hwid[, terms, drop = FALSE]
            lwr <- predictor + hwid
            upr <- predictor - hwid
        }
    }
    if (se.fit || interval != "none") {
        se <- sqrt(ip)
        if (type == "terms" && !is.null(terms) && !se.fit)
            se <- se[, terms, drop = FALSE]
    }
    if (missing(newdata) && !is.null(na.act <- object$na.action)) {
        predictor <- napredict(na.act, predictor)
        if (se.fit)
            se <- napredict(na.act, se)
    }
    if (type == "terms" && interval != "none") {
        if (missing(newdata) && !is.null(na.act)) {
            lwr <- napredict(na.act, lwr)
            upr <- napredict(na.act, upr)
        }
        list(fit = predictor, se.fit = se, lwr = lwr, upr = upr,
            df = df, residual.scale = sqrt(res.var))
    }
    else if (se.fit)
        list(fit = predictor, se.fit = se, df = df, residual.scale = sqrt(res.var))
    else predictor
}



getQr <- function(x, ...){
  if (is.null(r <- x$qr))
        stop("lm object does not have a proper 'qr' component.\n Rank zero or should not have used lm(.., qr=FALSE).")
  r
}
