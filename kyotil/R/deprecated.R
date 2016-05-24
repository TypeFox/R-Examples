getModeFromLocfit=function (fit) {
    tmp=preplot(fit)
    tmp$xev[[1]][rank (tmp$fit) == length(tmp$fit)]
}

setSeedAsInR=function() {
# choose a rng so that results can be compared with that from C program
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = NULL)
    set.seed(1, kind = "Marsaglia-Multicarry")
    myprint(runif(1))
    myprint("the number above should be 0.006153224")
}

#add CI to summary
mysummary <- function(object, ...) UseMethod("mysummary") 

#add prediction interval to summary
mypredict <- function(object, ...) UseMethod("mypredict") 

mysummary.lm = function (object, correlation = FALSE, symbolic.cor = FALSE, ...) 
{
    z <- object
    p <- z$rank
    if (p == 0) {
        r <- z$residuals
        n <- length(r)
        w <- z$weights
        if (is.null(w)) {
            rss <- sum(r^2)
        }
        else {
            rss <- sum(w * r^2)
            r <- sqrt(w) * r
        }
        resvar <- rss/(n - p)
        ans <- z[c("call", "terms")]
        class(ans) <- "summary.lm"
        ans$aliased <- is.na(coef(object))
        ans$residuals <- r
        ans$df <- c(0, n, length(ans$aliased))
        ans$coefficients <- matrix(NA, 0, 6)
        dimnames(ans$coefficients) <- list(NULL, c("Estimate", 
            "Std. Error", "lower CI", "higher CI", "t value", "Pr(>|t|)"))
        ans$sigma <- sqrt(resvar)
        ans$r.squared <- ans$adj.r.squared <- 0
        return(ans)
    }
    Qr <- object$qr
    if (is.null(z$terms) || is.null(Qr)) 
        stop("invalid 'lm' object:  no terms nor qr component")
    n <- NROW(Qr$qr)
    rdf <- n - p
    if (rdf != z$df.residual) 
        warning("inconsistent residual degrees of freedom. -- please report!")
    p1 <- 1:p
    r <- z$residuals
    f <- z$fitted
    w <- z$weights
    if (is.null(w)) {
        mss <- if (attr(z$terms, "intercept")) 
            sum((f - mean(f))^2)
        else sum(f^2)
        rss <- sum(r^2)
    }
    else {
        mss <- if (attr(z$terms, "intercept")) {
            m <- sum(w * f/sum(w))
            sum(w * (f - m)^2)
        }
        else sum(w * f^2)
        rss <- sum(w * r^2)
        r <- sqrt(w) * r
    }
    resvar <- rss/rdf
    R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
    se <- sqrt(diag(R) * resvar)
    est <- z$coefficients[Qr$pivot[p1]]
    #youyi
    ci.l = est - qt(.975, rdf, lower.tail=TRUE)*se
    ci.r = est + qt(.975, rdf, lower.tail=TRUE)*se
    tval <- est/se
    ans <- z[c("call", "terms")]
    ans$residuals <- r
    ans$coefficients <- cbind(est, se, tval, 2 * pt(abs(tval), 
        rdf, lower.tail = FALSE), ci.l, ci.r)
    dimnames(ans$coefficients) <- list(names(z$coefficients)[Qr$pivot[p1]], 
        c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "CI left", "CI right"))
    ans$aliased <- is.na(coef(object))
    ans$sigma <- sqrt(resvar)
    ans$df <- c(p, rdf, NCOL(Qr$qr))
    if (p != attr(z$terms, "intercept")) {
        df.int <- if (attr(z$terms, "intercept")) 
            1
        else 0
        ans$r.squared <- mss/(mss + rss)
        ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n - 
            df.int)/rdf)
        ans$fstatistic <- c(value = (mss/(p - df.int))/resvar, 
            numdf = p - df.int, dendf = rdf)
    }
    else ans$r.squared <- ans$adj.r.squared <- 0
    ans$cov.unscaled <- R
    dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1, 
        1)]
    if (correlation) {
        ans$correlation <- (R * resvar)/outer(se, se)
        dimnames(ans$correlation) <- dimnames(ans$cov.unscaled)
        ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- "summary.lm"
    ans
}

#mysummary.ols=function (x, digits = 4, long = FALSE, ...) 
#{
#    oldopt <- options(digits = digits)
#    on.exit(options(oldopt))
#    cat("\n")
#    cat("Linear Regression Model\n\n")
#    dput(x$call)
#    cat("\n")
#    if (!is.null(z <- x$na.action)) 
#        naprint(z)
#    stats <- x$stats
#    if (lst <- length(stats)) {
#        if (.R.) {
#            cstats <- character(lst)
#            names(cstats) <- names(stats)
#            for (i in 1:lst) cstats[i] <- format(stats[i])
#            print(cstats, quote = FALSE)
#        }
#        else print(x$stats)
#        cat("\n")
#    }
#    pen <- length(x$penalty.matrix) > 0
#    resid <- x$residuals
#    n <- length(resid)
#    p <- length(x$coef) - (names(x$coef)[1] == "Intercept")
#    if (length(x$stats) == 0) 
#        cat("n=", n, "   p=", p, "\n\n", sep = "")
#    ndf <- x$stats["d.f."]
#    df <- c(ndf, n - ndf - 1, ndf)
#    r2 <- x$stats["R2"]
#    sigma <- x$stats["Sigma"]
#    rdf <- df[2]
#    if (rdf > 5) {
#        cat("Residuals:\n")
#        if (length(dim(resid)) == 2) {
#            rq <- apply(t(resid), 1, quantile)
#            dimnames(rq) <- list(c("Min", "1Q", "Median", "3Q", 
#                "Max"), dimnames(resid)[[2]])
#        }
#        else {
#            rq <- quantile(resid)
#            names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
#        }
#        print(rq, digits = digits, ...)
#    }
#    else if (rdf > 0) {
#        cat("Residuals:\n")
#        print(resid, digits = digits, ...)
#    }
#    if (nsingular <- df[3] - df[1]) 
#        cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n", 
#            sep = "")
#    else cat("\nCoefficients:\n")
#    se <- sqrt(diag(x$var))
#    z <- x$coefficients/se
#    P <- 2 * (1 - pt(abs(z), rdf))
#    ci.l = x$coefficients - qt(.975, rdf, lower.tail=TRUE)*se
#    ci.r = x$coefficients + qt(.975, rdf, lower.tail=TRUE)*se
#    co <- cbind(x$coefficients, se, z, P, ci.l, ci.r)
#    dimnames(co) <- list(names(x$coefficients), c("Value", "Std. Error", 
#        "t", "Pr(>|t|)", "CI left", "CI right"))
#    print(co)
#    if (pen) 
#        cat("\n")
#    else cat("\nResidual standard error:", format(signif(sigma, 
#        digits)), "on", rdf, "degrees of freedom\n")
#    rsqa <- 1 - (1 - r2) * (n - 1)/rdf
#    if (length(x$stats) == 0) 
#        cat("Multiple R-Squared:", format(signif(r2, digits)), 
#            " ")
#    cat("Adjusted R-Squared:", format(signif(rsqa, digits)), 
#        "\n")
#    if (!pen) {
#        if (long && p > 0) {
#            correl <- diag(1/se) %*% x$var %*% diag(1/se)
#            dimnames(correl) <- dimnames(x$var)
#            cat("\nCorrelation of Coefficients:\n")
#            ll <- lower.tri(correl)
#            correl[ll] <- format(round(correl[ll], digits), ...)
#            correl[!ll] <- ""
#            print(correl[-1, -(p + 1), drop = FALSE], quote = FALSE, 
#                digits = digits, ...)
#        }
#    }
#    cat("\n")
#    invisible()
#}
#
#summary.ols = function (ols1) {
#    print(ols1)
#}

mysummary.glm = function (object, dispersion = NULL, correlation = FALSE, symbolic.cor = FALSE, ...) 
{
    est.disp <- FALSE
    df.r <- object$df.residual
    if (is.null(dispersion)) 
        dispersion <- if (any(object$family$family == c("poisson", 
            "binomial"))) 
            1
        else if (df.r > 0) {
            est.disp <- TRUE
            if (any(object$weights == 0)) 
                warning("observations with zero weight ", "not used for calculating dispersion")
            sum(object$weights * object$residuals^2)/df.r
        }
        else Inf
    p <- object$rank
    if (p > 0) {
        p1 <- 1:p
        Qr <- object$qr
        aliased <- is.na(coef(object))
        coef.p <- object$coefficients[Qr$pivot[p1]]
        covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
        dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
        covmat <- dispersion * covmat.unscaled
        var.cf   <- diag(covmat)
        s.err <- sqrt(var.cf)
        tvalue <- coef.p/s.err
        
        #youyi
        ci.l = coef.p - qt(.975, df.r, lower.tail=TRUE)*s.err
        ci.r = coef.p + qt(.975, df.r, lower.tail=TRUE)*s.err
    
        dn <- c("Estimate", "Std. Error")
        if (!est.disp) {
            pvalue <- 2 * pnorm(-abs(tvalue))
            coef.table <- cbind(coef.p, s.err, tvalue, pvalue, ci.l, ci.r)
            dimnames(coef.table) <- list(names(coef.p), c(dn, 
                "z value", "Pr(>|z|)", "lower CI", "upper CI"))
        }
        else if (df.r > 0) {
            pvalue <- 2 * pt(-abs(tvalue), df.r)
            coef.table <- cbind(coef.p, s.err, tvalue, pvalue, ci.l, ci.r)
            dimnames(coef.table) <- list(names(coef.p), c(dn, 
                "t value", "Pr(>|t|)", "lower CI", "upper CI"))
        }
        else {
            coef.table <- cbind(coef.p, Inf)
            dimnames(coef.table) <- list(names(coef.p), dn)
        }
        df.f <- NCOL(Qr$qr)
    }
    else {
        coef.table <- matrix(, 0, 4)
        dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error", 
            "t value", "Pr(>|t|)"))
        covmat.unscaled <- covmat <- matrix(, 0, 0)
        aliased <- is.na(coef(object))
        df.f <- length(aliased)
    }
    ans <- c(object[c("call", "terms", "family", "deviance", 
        "aic", "contrasts", "df.residual", "null.deviance", "df.null", 
        "iter")], list(deviance.resid = residuals(object, type = "deviance"), 
        coefficients = coef.table, aliased = aliased, dispersion = dispersion, 
        df = c(object$rank, df.r, df.f), cov.unscaled = covmat.unscaled, 
        cov.scaled = covmat))
    if (correlation && p > 0) {
        dd <- sqrt(diag(covmat.unscaled))
        ans$correlation <- covmat.unscaled/outer(dd, dd)
        ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- "summary.glm"
    return(ans)
}

mypredict.glm = function (object, newdata = NULL, type = c("link", "response", 
    "terms"), se.fit = FALSE, dispersion = NULL, terms = NULL, 
    na.action = na.pass, ...) 
{
    type <- match.arg(type)
    na.act <- object$na.action
    object$na.action <- NULL
    if (!se.fit) {
        if (missing(newdata)) {
            pred <- switch(type, link = object$linear.predictors, 
                response = object$fitted, terms = predict.lm(object, 
                  se.fit = se.fit, scale = 1, type = "terms", 
                  terms = terms))
            if (!is.null(na.act)) 
                pred <- napredict(na.act, pred)
        }
        else {
            pred <- predict.lm(object, newdata, se.fit, scale = 1, 
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
        pred <- predict.lm(object, newdata, se.fit, scale = residual.scale, 
            type = ifelse(type == "link", "response", type), 
            terms = terms, na.action = na.action, interval="confidence")#youyi adds interval here
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

mypredict.lm=function (object, newdata, se.fit = FALSE, scale = NULL, df = Inf, 
    interval = c("none", "confidence", "prediction"), level = 0.95, 
    type = c("response", "terms"), terms = NULL, na.action = na.pass, 
    ...) 
{
    tt <- terms(object)
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
        X <- model.matrix(Terms, m, contrasts = object$contrasts)
        offset <- if (!is.null(off.num <- attr(tt, "offset"))) 
            eval(attr(tt, "variables")[[off.num + 1]], newdata)
        else if (!is.null(object$offset)) 
            eval(object$call$offset, newdata)
        mmDone <- FALSE
    }
    n <- length(object$residuals)
    p <- object$rank
    p1 <- seq(len = p)
    piv <- object$qr$pivot[p1]
    if (p < ncol(X) && !(missing(newdata) || is.null(newdata))) 
        warning("prediction from a rank-deficient fit may be misleading")
    beta <- object$coefficients
    predictor <- drop(X[, piv, drop = FALSE] %*% beta[piv])
    if (!is.null(offset)) 
        predictor <- predictor + offset
    interval <- match.arg(interval)
    type <- match.arg(type)
    if (se.fit || interval != "none") {
        res.var <- if (is.null(scale)) {
            r <- object$residuals
            w <- object$weights
            rss <- sum(if (is.null(w)) 
                r^2
            else r^2 * w)
            df <- n - p
            rss/df
        }
        else scale^2
        if (type != "terms") {
            if (p > 0) {
                XRinv <- if (missing(newdata) && is.null(w)) 
                  qr.Q(object$qr)[, p1, drop = FALSE]
                else X[, piv] %*% qr.solve(qr.R(object$qr)[p1, 
                  p1])
                ip <- drop(XRinv^2 %*% rep(res.var, p))
            }
            else ip <- rep(0, n)
        }
    }
    if (type == "terms") {
        if (!mmDone) {
            mm <- model.matrix(object)
            mmDone <- TRUE
        }
        aa <- attr(mm, "assign")
        ll <- attr(tt, "term.labels")
        if (attr(tt, "intercept") > 0) 
            ll <- c("(Intercept)", ll)
        aaa <- factor(aa, labels = ll)
        asgn <- split(order(aa), aaa)
        hasintercept <- attr(tt, "intercept") > 0
        if (hasintercept) {
            asgn$"(Intercept)" <- NULL
            if (!mmDone) {
                mm <- model.matrix(object)
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
                Rinv <- qr.solve(qr.R(object$qr)[p1, p1])
            }
            if (hasintercept) 
                X <- sweep(X, 2, avx)
            unpiv <- rep.int(0, NCOL(X))
            unpiv[piv] <- p1
            for (i in seq(1, nterms, length = nterms)) {
                iipiv <- asgn[[i]]
                ii <- unpiv[iipiv]
                iipiv[ii == 0] <- 0
                predictor[, i] <- if (any(iipiv) > 0) 
                  X[, iipiv, drop = FALSE] %*% beta[iipiv]
                else 0
                if (se.fit || interval != "none") 
                  ip[, i] <- if (any(iipiv) > 0) 
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
            predictor <- ip <- matrix(0, n, 0)
        }
        attr(predictor, "constant") <- if (hasintercept) 
            termsconst
        else 0
    }
    if (interval != "none") {
        tfrac <- qt((1 - level)/2, df)
        hwid <- tfrac * switch(interval, confidence = sqrt(ip), 
            prediction = sqrt(ip + res.var))
        if (type != "terms") {
            predictor <- cbind(predictor, predictor + hwid %o% 
                c(1, -1))
            colnames(predictor) <- c("fit", "lwr", "upr")
        }
        else {
            lwr <- predictor + hwid
            upr <- predictor - hwid
        }
    }
    if (se.fit || interval != "none") 
        se <- sqrt(ip)
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

# predict.gam returns a list of fit and se.fit. This function transforms it to a matrix: 1st col is fit, 2nd col is LI, 3rd col is UI, 4th col is se.fit
mypredict.gam = function (...) {
    pred = predict(...)
    x=matrix (0, length (pred$fit), 4)
    x[,1]=pred$fit
    x[,2]=pred$fit - 2 * pred$se.fit
    x[,3]=pred$fit + 2 * pred$se.fit
    x[,4]=pred$se.fit
    dimnames (x) = list (NULL, list("prediction", "LI", "UI", "se"))
    x
}
