#  Modification of summary.glm from the stats package for R.
#
#  Copyright (C) 1995-2005 The R Core Team
#  Copyright (C) 2005, 2006, 2010, 2015 Heather Turner
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

summary.gnm <- function (object, dispersion = NULL, correlation = FALSE,
                         symbolic.cor = FALSE, with.eliminate = FALSE, ...)
{
    est.disp <- (!object$family$family %in% c("poisson", "binomial") &&
                 is.null(dispersion) && object$df.residual > 0)
    coefs <- parameters(object)
    if (with.eliminate) coefs <- c(attr(coef(object), "eliminated"), coefs)
    if (object$rank > 0) {
        cov.scaled <- vcov(object, dispersion = dispersion,
                           with.eliminate = with.eliminate)
        ## non-eliminated par only
        if (nrow(cov.scaled)) {
            estimable <- checkEstimable(object, ...)
            estimable[is.na(estimable)] <- FALSE
        }
        if (is.matrix(cov.scaled))
            sterr <- sqrt(diag(cov.scaled))
        else
            sterr <- diag(cov.scaled)
        if (length(sterr)) is.na(sterr[!estimable]) <- TRUE
        if (with.eliminate){
            ## check estimability of eliminated coefficients
            X <- cbind(1, model.matrix(object)[,!is.na(coef(object))])
            estimable2 <- sapply(split(1:nrow(X), object$eliminate),
                                 function(i) {
                                     quickRank(X[i, , drop = FALSE]) ==
                                         quickRank(X[i, -1, drop = FALSE]) + 1})
            sterr <- c(ifelse(estimable2,
                              sqrt(attr(cov.scaled, "varElim")), NA),
                       sterr)
        }
        tvalue <- coefs/sterr
        dn <- c("Estimate", "Std. Error")
        if (!est.disp) {
            pvalue <- 2 * pnorm(-abs(tvalue))
            coef.table <- cbind(coefs, sterr, tvalue, pvalue)
            dimnames(coef.table) <- list(names(coefs),
                                         c(dn, "z value", "Pr(>|z|)"))
        }
        else if (object$df.residual > 0) {
            pvalue <- 2 * pt(-abs(tvalue), object$df.residual)
            coef.table <- cbind(coefs, sterr, tvalue, pvalue)
            dimnames(coef.table) <- list(names(coefs),
                                         c(dn, "t value", "Pr(>|t|)"))
        }
        else {
            coef.table <- cbind(coefs, Inf)
            dimnames(coef.table) <- list(names(coefs), dn)
        }
    }
    else {
        coef.table <- matrix(, 0, 4)
        dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error",
            "t value", "Pr(>|t|)"))
        cov.scaled <- matrix(, 0, 0)
    }
    df.f <- nrow(coef.table)
    non.elim <- seq(object$coef) + nlevels(object$eliminate)*with.eliminate
    elim <- seq(length.out = nlevels(object$eliminate)*with.eliminate)
    ans <- c(object[c("call", "ofInterest", "family", "deviance", "aic",
                      "df.residual", "iter")],
             list(deviance.resid = residuals(object, type = "deviance"),
                  coefficients = coef.table[non.elim, , drop = FALSE],
                  eliminated = coef.table[elim, , drop = FALSE],
                  dispersion = attr(cov.scaled, "dispersion"),
                  df = c(object$rank, object$df.residual, df.f),
                  cov.scaled = as.matrix(cov.scaled)))
    if (correlation & object$rank > 0) {
        dd <- sqrt(diag(cov.scaled))
        ans$correlation <- cov.scaled/outer(dd, dd)
        ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- "summary.gnm"
    ans
}
