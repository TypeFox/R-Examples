#  This fits a glm with eliminated factor, and should be much quicker
#  than glm.fit when the number of levels of the eliminated factor is large.
#
#  Copyright (C) 2009, 2010, 2012 David Firth and Heather Turner
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

glm.fit.e <- function(
    x, y,
    weights = rep(1, NROW(y)),
    start = NULL,
    etastart = NULL,
    mustart = NULL,
    offset = rep(0, NROW(y)),
    family = gaussian(),
    control = glm.control(), ## only for compatibility with glm.fit
    intercept = TRUE, ## only for compatibility with glm.fit
    eliminate = NULL,  ## alternatively a factor
    ridge = 1e-8,
    coefonly = FALSE)
{
    if (is.null(eliminate)) { ## just revert to glm.fit
        ## can make a difference in timing!
        tmp <- glm.fit(x, y, weights = weights, start = start,
                           etastart = etastart, mustart = mustart,
                           offset = offset, family = family,
                           control = control,
                           intercept = intercept)
        if (coefonly) return(tmp$coef)
        else return(tmp)
    }
    ##  The rest handles the case of an eliminated factor
    names(y) <- rownames(x) <- NULL
    nobs <- NROW(y)
    non.elim <- ncol(x)
    if (is.null(weights))
        weights <- rep.int(1, nobs)
    if (is.null(offset))
        offset <- rep.int(0, nobs)
    link <- family$linkfun
    linkinv <- family$linkinv
    linkder <- family$mu.eta
    variance <- family$variance
    dev.resids <- family$dev.resids
    aic <- family$aic
    ## sort data to help compute group means quickly
    ord <- order(xtfrm(eliminate))
    if (ordTRUE <- !identical(ord, seq(eliminate))) {
        y <- as.numeric(y[ord])
        weights <- weights[ord]
        offset <- offset[ord]
        if (non.elim) x <- x[ord, , drop = FALSE]
        eliminate <- eliminate[ord]
    }
    size <- tabulate(eliminate)
    end <- cumsum(size)
    nelim <- rank <- nlevels(eliminate)
    elim <- seq.int(nelim)
    if (is.null(start)) { # use either y or etastart or mustart
        if (!is.null(etastart)) mustart <- linkinv(etastart)
        if (!is.null(mustart)) z <- mustart
        else z <- y
        elim.means <- grp.sum(z, end)/size
        os.by.level <- link(0.999 * elim.means + 0.001 * mean(z)) -
            grp.sum(offset, end)/size
    } else os.by.level <- start[elim]
    os.vec <- os.by.level[eliminate]
    eta.stored <- eta <- offset + os.vec
    mu <- linkinv(eta)
    mu.eta <- linkder(eta)
    z <- eta - offset + (y - mu) / mu.eta
    w <- weights * (mu.eta)^2/variance(mu)
    counter <- 0
    devold <- 0
    if (intercept) x <- x[, -1, drop = FALSE] #non-null eliminate
    if (non.elim) {
        ## sweeps needed to get the rank right
        subtracted <- rowsum.default(x, eliminate, reorder = FALSE)/size
        x <- x - subtracted[eliminate, , drop = FALSE]
        ## initial fit to drop aliased columns
        model <- lm.wfit(x, z, w, offset = os.vec)
        full.theta <- model$coefficients
        eta <- model$fitted + offset
        rank <- model$rank + nelim
        rm(model)
        mu <- linkinv(eta)
        mu.eta <- linkder(eta)
        z <- eta - offset + (y - mu) / mu.eta
        w <- weights * (mu.eta)^2/variance(mu)
        est <- !is.na(full.theta)
        x <- x[, est, drop = FALSE]
        theta <- full.theta[est]
    }
    Z <- cbind(z, x)
    I1 <- numeric(ncol(Z))
    I1[1] <- 1
    for (i in 1:control$maxit) {
        ## try without scaling etc - already of full rank
        Tvec <- sqrt(grp.sum(w, end))
        Umat <- rowsum.default(w * Z, eliminate, reorder = FALSE)
        Umat <- Umat/Tvec
        Wmat <- crossprod(sqrt(w) * Z)
        diag(Wmat) <- diag(Wmat) + ridge
        Qi <- solve(Wmat - crossprod(Umat), I1)
        theta <- -Qi[-1]/Qi[1]
        os.by.level <- ((Umat %*% Qi)/Qi[1])/Tvec
        if (non.elim) eta <- drop(x %*% theta + offset + os.by.level[eliminate])
        else eta <- offset + os.by.level[eliminate]
        mu <- linkinv(eta)
        dev <- sum(dev.resids(y, mu, weights))
        if (control$trace)
            cat("Deviance =", dev, "Iterations -", i,
                "\n")
        if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
            conv <- TRUE
            break
        }
        devold <- dev
        mu.eta <- linkder(eta)
        Z[,1] <- eta - offset + (y - mu) / mu.eta
        w <- weights * (mu.eta)^2/variance(mu)
    }
    converged <- !(i == control$maxit)
    if (!converged) warning(paste("The convergence criterion was not met after",
                                  control$maxit, "iterations."))
    names(os.by.level) <- paste("(eliminate)", elim, sep = "")
    if (non.elim) {
        full.theta[est] <- theta
        os.by.level <- os.by.level - subtracted %*% naToZero(full.theta)
    }
    else full.theta <- numeric(0)
    if (ordTRUE) {
        reorder <- order(ord)
        y <- y[reorder]
        mu <- mu[reorder]
        eta <- eta[reorder]
        weights <- weights[reorder]
    }
    mu.eta <- linkder(eta)
    w <- weights * (mu.eta)^2/variance(mu)
    if (coefonly) return(structure(full.theta, eliminated = c(os.by.level)))
    aic.model <- aic(y, sum(weights > 0), mu, weights, dev) + 2 * rank
    eliminated <- structure(c(os.by.level), names = levels(eliminate))
    list(coefficients = structure(full.theta, eliminated = eliminated),
         residuals = (y - mu) / linkder(eta),
         fitted.values = mu,
         rank = rank,
         family = family,
         linear.predictors = eta,
         deviance = dev,
         aic = aic.model,
         iter = i,
         weights = w,
         prior.weights = weights,
         df.residual = nobs - sum(weights == 0) - rank,
         y = y,
         converged = converged)
    ##  NB: some components of the result of glm.fit are missing from this list
}


