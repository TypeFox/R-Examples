#  Modification of drop1.glm from the stats package for R.
#
#  Copyright (C) 1994-8 W. N. Venables and B. D. Ripley
#  Copyright (C) 1998-2005 The R Core Team
#  Copyright (C) 2005, 2010, 2013 Heather Turner
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

drop1.gnm <- function (object, scope, scale = 0, test = c("none", "Chisq",
    "F"), k = 2, ...)
{
    if (any(attr(terms(object), "type") != "Linear"))
        stop("add1 is not implemented for gnm objects with nonlinear terms.")
    x <- model.matrix(object)
    n <- nrow(x)
    asgn <- attr(x, "assign")
    tl <- attr(object$terms, "term.labels")
    if (missing(scope))
        scope <- drop.scope(object)
    else {
        if (!is.character(scope))
            scope <- attr(terms(update.formula(object, scope)),
                "term.labels")
        if (!all(match(scope, tl, 0L) > 0L))
            stop("scope is not a subset of term labels")
    }
    ndrop <- match(scope, tl)
    ns <- length(scope)
    rdf <- object$df.residual
    chisq <- object$deviance
    dfs <- numeric(ns)
    dev <- numeric(ns)
    y <- object$y
    if (is.null(y)) {
        y <- model.response(model.frame(object))
        if (!is.factor(y))
            storage.mode(y) <- "double"
    }
    wt <- object$prior.weights
    if (is.null(wt))
        wt <- rep.int(1, n)
    for (i in 1L:ns) {
        ii <- seq_along(asgn)[asgn == ndrop[i]]
        jj <- setdiff(seq(ncol(x)), ii)
        z <- glm.fit.e(x[, jj, drop = FALSE], y, wt, offset = object$offset,
            family = object$family, eliminate = object$eliminate)
        dfs[i] <- z$rank
        dev[i] <- z$deviance
    }
    scope <- c("<none>", scope)
    dfs <- c(object$rank, dfs)
    dev <- c(chisq, dev)
    dispersion <- if (is.null(scale) || scale == 0)
        summary(object, dispersion = NULL)$dispersion
    else scale
    fam <- object$family$family
    loglik <- if (fam == "gaussian") {
        if (scale > 0)
            dev/scale - n
        else n * log(dev/n)
    }
    else dev/dispersion
    aic <- loglik + k * dfs
    dfs <- dfs[1L] - dfs
    dfs[1L] <- NA
    aic <- aic + (extractAIC(object, k = k)[2L] - aic[1L])
    aod <- data.frame(Df = dfs, Deviance = dev, AIC = aic, row.names = scope,
        check.names = FALSE)
    if (all(is.na(aic)))
        aod <- aod[, -3]
    test <- match.arg(test)
    if (test == "Chisq") {
        dev <- pmax(0, loglik - loglik[1L])
        dev[1L] <- NA
        nas <- !is.na(dev)
        LRT <- if (dispersion == 1)
            "LRT"
        else "scaled dev."
        aod[, LRT] <- dev
        dev[nas] <- pchisq(dev[nas], aod$Df[nas], lower.tail = FALSE)
        aod[, "Pr(Chi)"] <- dev
    }
    else if (test == "F") {
        if (fam == "binomial" || fam == "poisson")
            warning(gettextf("F test assumes 'quasi%s' family",
                fam), domain = NA)
        dev <- aod$Deviance
        rms <- dev[1L]/rdf
        dev <- pmax(0, dev - dev[1L])
        dfs <- aod$Df
        rdf <- object$df.residual
        Fs <- (dev/dfs)/rms
        Fs[dfs < 1e-04] <- NA
        P <- Fs
        nas <- !is.na(Fs)
        P[nas] <- pf(Fs[nas], dfs[nas], rdf, lower.tail = FALSE)
        aod[, c("F value", "Pr(F)")] <- list(Fs, P)
    }
    head <- c("Single term deletions", "\nModel:", deparse(as.vector(formula(object))),
        if (!is.null(scale) && scale > 0) paste("\nscale: ",
            format(scale), "\n"))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}
