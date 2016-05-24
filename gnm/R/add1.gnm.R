#  Modification of add1.glm from the stats package for R.
#
#  Copyright (C) 1994-8 W. N. Venables and B. D. Ripley
#  Copyright (C) 1998-2005 The R Core Team
#  Copyright (C) 2005, 2010 Heather Turner
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

add1.gnm <- function (object, scope, scale = 0,
                      test = c("none", "Chisq", "F"), x = NULL, k = 2, ...)
{
    if (any(attr(terms(object), "type") != "Linear"))
        stop("add1 is not implemented for gnm objects with nonlinear terms.")
    Fstat <- function(table, rdf) {
        dev <- table$Deviance
        df <- table$Df
        diff <- pmax(0, (dev[1L] - dev)/df)
        Fs <- (diff/df)/(dev/(rdf - df))
        Fs[df < .Machine$double.eps] <- NA
        P <- Fs
        nnas <- !is.na(Fs)
        P[nnas] <- pf(Fs[nnas], df[nnas], rdf - df[nnas],
            lower.tail = FALSE)
        list(Fs = Fs, P = P)
    }
    if (!is.character(scope))
        scope <- add.scope(object, update.formula(object, scope))
    if (!length(scope))
        stop("no terms in scope for adding to object")
    oTerms <- attr(object$terms, "term.labels")
    int <- attr(object$terms, "intercept")
    ns <- length(scope)
    dfs <- dev <- numeric(ns + 1)
    names(dfs) <- names(dev) <- c("<none>", scope)
    add.rhs <- paste(scope, collapse = "+")
    add.rhs <- eval(parse(text = paste("~ . +", add.rhs)))
    new.form <- update.formula(object, add.rhs)
    Terms <- terms(new.form)
    y <- object$y
    if (is.null(x)) {
        fc <- object$call
        fc$formula <- Terms
        fob <- list(call = fc, terms = Terms)
        class(fob) <- oldClass(object)
        m <- model.frame(fob, xlev = object$xlevels)
        offset <- model.offset(m)
        wt <- model.weights(m)
        x <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
        oldn <- length(y)
        y <- model.response(m)
        if (!is.factor(y))
            storage.mode(y) <- "double"
        if (NCOL(y) == 2) {
            n <- y[, 1] + y[, 2]
            y <- ifelse(n == 0, 0, y[, 1]/n)
            if (is.null(wt))
                wt <- rep.int(1, length(y))
            wt <- wt * n
        }
        newn <- length(y)
        if (newn < oldn)
            warning(gettextf("using the %d/%d rows from a combined fit",
                newn, oldn), domain = NA)
    }
    else {
        wt <- object$prior.weights
        offset <- object$offset
    }
    n <- nrow(x)
    if (is.null(wt))
        wt <- rep.int(1, n)
    Terms <- attr(Terms, "term.labels")
    asgn <- attr(x, "assign")
    ousex <- match(asgn, match(oTerms, Terms), 0L) > 0L
    if (int)
        ousex[1L] <- TRUE
    X <- x[, ousex, drop = FALSE]
    z <- glm.fit.e(X, y, wt, offset = offset, family = object$family,
                   eliminate = object$eliminate)
    dfs[1L] <- z$rank
    dev[1L] <- z$deviance
    sTerms <- sapply(strsplit(Terms, ":", fixed = TRUE), function(x) paste(sort(x),
        collapse = ":"))
    for (tt in scope) {
        stt <- paste(sort(strsplit(tt, ":")[[1L]]), collapse = ":")
        usex <- match(asgn, match(stt, sTerms), 0L) > 0L
        X <- x[, usex | ousex, drop = FALSE]
        z <- glm.fit.e(X, y, wt, offset = offset, family = object$family,
                       eliminate = object$eliminate)
        dfs[tt] <- z$rank
        dev[tt] <- z$deviance
    }
    if (scale == 0)
        dispersion <- summary(object, dispersion = NULL)$dispersion
    else dispersion <- scale
    fam <- object$family$family
    if (fam == "gaussian") {
        if (scale > 0)
            loglik <- dev/scale - n
        else loglik <- n * log(dev/n)
    }
    else loglik <- dev/dispersion
    aic <- loglik + k * dfs
    aic <- aic + (extractAIC(object, k = k)[2L] - aic[1L])
    dfs <- dfs - dfs[1L]
    dfs[1L] <- NA
    aod <- data.frame(Df = dfs, Deviance = dev, AIC = aic, row.names = names(dfs),
        check.names = FALSE)
    if (all(is.na(aic)))
        aod <- aod[, -3]
    test <- match.arg(test)
    if (test == "Chisq") {
        dev <- pmax(0, loglik[1L] - loglik)
        dev[1L] <- NA
        LRT <- if (dispersion == 1)
            "LRT"
        else "scaled dev."
        aod[, LRT] <- dev
        nas <- !is.na(dev)
        dev[nas] <- pchisq(dev[nas], aod$Df[nas], lower.tail = FALSE)
        aod[, "Pr(Chi)"] <- dev
    }
    else if (test == "F") {
        if (fam == "binomial" || fam == "poisson")
            warning(gettextf("F test assumes quasi%s family",
                fam), domain = NA)
        rdf <- object$df.residual
        aod[, c("F value", "Pr(F)")] <- Fstat(aod, rdf)
    }
    head <- c("Single term additions", "\nModel:", deparse(as.vector(formula(object))),
        if (scale > 0) paste("\nscale: ", format(scale), "\n"))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}
