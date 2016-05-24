#  Copyright (C) 2012 Yohan Chalabi
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License as
#  published by the Free Software Foundation; either version 2 or 3 of
#  the License (at your option).
#
#  This program is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

test.dgl <- function() {

    # since dgl is based on pgl and qdgl and that these functions have
    # been checked, we can simply implementation Rdgl with them:
    Rdgl <- function(x, pars) {
        p <- pgl(x, pars)
        1/qdgl(p, pars)
    }

    tol <- 5 * .Machine$double.eps

    ###################################
    # generate different sets of parameter values
    iqr <- c(1, 2)
    med <- c(-1, 0, 1)
    chi <- c(-.75, -.5, -.25, 0, .25, .5, .75)    # take care to include
    xi <- unique(sort(.5 * (1. + c(chi, -chi))))  # special cases

    params <- expand.grid(med, iqr, chi, xi)

    params <- rbind(params,
                    c(0, 1, 0, .5 - 1/sqrt(5)), # special case
                    c(0, 1, 0, .5 - 2/sqrt(17)),
                    c(0, 1, 0, .5),
                    c(0, 1,  1, 0),
                    c(0, 1, -1, 0))

    # special case with different location and scale
    params <- rbind(params,
                    c(-1, 2, 0, .5 - 1/sqrt(5)),
                    c(-1, 2, 0, .5 - 2/sqrt(17)),
                    c(-1, 2, 0, .5),
                    c(-1, 2,  1, 0),
                    c(-1, 2, -1, 0))
    params <- rbind(params,
                    c(1, 2, 0, .5 - 1/sqrt(5)),
                    c(1, 2, 0, .5 - 2/sqrt(17)),
                    c(1, 2, 0, .5),
                    c(1, 2,  1, 0),
                    c(1, 2, -1, 0))

    ###################################
    # check if density integrates to one
    for (i in seq.int(nrow(params))) {
        pars <- as.numeric(params[i, ])
        range <- qgl(c(0, 1), pars)
        den <- integrate(dgl, range[1], range[2], med = pars,
                         subdivisions = 1e4)
        checkTrue(den$value && (den$abs.error < 1e-3))
    }

    ###################################
    # check if we have the same result as with R implementation
    for (i in seq.int(nrow(params))) {
        pars <- as.numeric(params[i, ])
        p <- seq(0, 1, len = 1e3)
        x <- qgl(p, pars)
        checkEquals(dgl(x, pars), Rdgl(x, pars))
    }

    ###################################
    # test included distribution

    p <- seq(0, 1, len = 1e3)

    # uniform
    a <- 1; b <- 5
    pars <- c(.5 * (a + b), .5 * (b - a), 0, .5 - 1/sqrt(5))

    x <- qunif(p, min = a, max = b)
    d <- dunif(x, min = a, max = b)
    gd <- dgl(x, pars)
    checkEquals(d, gd, tol)

    pars <- c(.5 * (a + b), .5 * (b - a), 0, .5 - 2/sqrt(17))
    gd <- dgl(x, pars)
    checkEquals(d, gd, tol)

    # logistic
    loc <- 2; scale <- 3
    q <- qlogis(p, loc, scale)
    med <- loc; iqr <- scale * log(9); chi <- 0; xi <- .5
    checkEquals(q, qgl(p, med, iqr, chi, xi), tol)

    # exponential
    rate <- 3
    q <- qexp(p, rate)
    med <- log(2)/rate; iqr <- log(3)/rate; chi <- 1; xi <- 0
    checkEquals(q, qgl(p, med, iqr, chi, xi), tol)
    med <- -log(2)/rate; iqr <- log(3)/rate; chi <- -1; xi <- 0
    checkEquals(q, -qgl(1-p, med, iqr, chi, xi), tol)

    ###################################
    # test with NaNs in vector
    for (i in 1:nrow(params)) {
        pars <- as.numeric(params[i, ])
        p <- c(0, NaN, .25, NaN, .5, NaN, .75, NaN, 1, NaN)
        x <- qgl(p, pars)
        d <- dgl(x, pars)
        # check if we really have NaNs back with pgl
        checkIdentical(x[c(FALSE, TRUE)], rep(NaN, 5))
        # and check that other values were not changed because of the
        # NaNs
        idx <- c(TRUE, FALSE)
        checkEquals(dgl(x[idx], pars), d[idx], tol = tol)
    }

    ###################################
    # test with NaNs in vector
    for (i in 1:nrow(params)) {
        pars <- as.numeric(params[i, ])
        p <- c(0, NA, .25, NA, .5, NA, .75, NA, 1, NA)
        x <- qgl(p, pars)
        d <- dgl(x, pars)
        # check if we really have NAs back with pgl
        checkIdentical(d[c(FALSE, TRUE)], rep(NA_real_, 5))
        # and check that other values were not changed because of the
        # NAs
        idx <- c(TRUE, FALSE)
        checkEquals(dgl(x[idx], pars), d[idx], tol = tol)
    }

    ###################################
    # test with values larger  or smaller than distribution range
    for (i in 1:nrow(params)) {

        pars <- as.numeric(params[i, ])
        xmin <- qgl(0, pars)
        xmax <- qgl(1, pars)

        if (is.finite(xmin)) {
            x <- c(if (xmin < 0) xmin * 1.2 else xmin * .8, xmin)
            d <- dgl(x, pars)
            checkEquals(d, c(0, 1/qdgl(0, pars)), tol)
        }

        if (is.finite(xmax)) {
            x <- c(xmax, if (xmax < 0) xmax * .8 else xmax * 1.2)
            d <- dgl(x, pars)
            checkEquals(d, c(1/qdgl(1, pars), 0), tol)
        }
    }

    ###################################
    # test invalid options that should produce NaNs and warnings

    oo <- options(warn = -1)
    on.exit(options(oo))

    checkTrue(all(is.nan(dgl(c(0, 0, 0), iqr = 0))))
    checkTrue(all(is.nan(dgl(c(0, 0, 0), iqr = -1))))
    checkTrue(all(is.nan(dgl(c(0, 0, 0), chi = -1))))
    checkTrue(all(is.nan(dgl(c(0, 0, 0), chi = 1))))
    checkTrue(all(is.nan(dgl(c(0, 0, 0), chi = 2))))
    checkTrue(all(is.nan(dgl(c(0, 0, 0), chi = -2))))
    checkTrue(all(is.nan(dgl(c(0, 0, 0), xi = 0))))
    checkTrue(all(is.nan(dgl(c(0, 0, 0), xi = 1))))
    checkTrue(all(is.nan(dgl(c(0, 0, 0), xi = -1))))
    checkTrue(all(is.nan(dgl(c(0, 0, 0), xi = 2))))

    options(warn = 2)

    checkException(dgl(c(0, 0, 0), iqr = 0), silent = TRUE)
    checkException(dgl(c(0, 0, 0), iqr = -1), silent = TRUE)
    checkException(dgl(c(0, 0, 0), chi = -1), silent = TRUE)
    checkException(dgl(c(0, 0, 0), chi = 1), silent = TRUE)
    checkException(dgl(c(0, 0, 0), chi = 2), silent = TRUE)
    checkException(dgl(c(0, 0, 0), chi = -2), silent = TRUE)
    checkException(dgl(c(0, 0, 0), xi = 0), silent = TRUE)
    checkException(dgl(c(0, 0, 0), xi = 1), silent = TRUE)
    checkException(dgl(c(0, 0, 0), xi = -1), silent = TRUE)
    checkException(dgl(c(0, 0, 0), xi = 2), silent = TRUE)

    options(oo)

}
