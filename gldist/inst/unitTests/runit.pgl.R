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

test.pgl <- function() {

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
    # test precision with a large vector of probabilities
    p <- seq(0, 1, len = 1e3)
    testP <- apply(params, 1, function(pars) {
        q <- qgl(p, pars)
        pp <- pgl(q, pars)
        max(abs((p - pp))) <= tol
    })
    checkTrue(all(testP))

    ###################################
    # test with NaNs in vector of probabilities
    for (i in 1:nrow(params)) {
        pars <- as.numeric(params[i, ])
        p <- c(0, NaN, .25, NaN, .5, NaN, .75, NaN, 1, NaN)
        q <- qgl(p, pars)
        # check if we really have NaNs back with pgl
        checkIdentical(q[c(FALSE, TRUE)], rep(NaN, 5))
        # and check that other values were not changed because of the
        # NaNs
        idx <- c(TRUE, FALSE)
        checkEquals(pgl(q[idx], pars), p[idx], tol = tol)
    }

    ###################################
    # test with NaNs in vector of probabilities
    for (i in 1:nrow(params)) {
        pars <- as.numeric(params[i, ])
        p <- c(0, NA, .25, NA, .5, NA, .75, NA, 1, NA)
        q <- qgl(p, pars)
        # check if we really have NAs back with pgl
        checkIdentical(q[c(FALSE, TRUE)], rep(NA_real_, 5))
        # and check that other values were not changed because of the
        # NAs
        idx <- c(TRUE, FALSE)
        checkEquals(pgl(q[idx], pars), p[idx], tol = tol)
    }

    ###################################
    # test values that are equal or smaller than 0 quantile
    for (i in seq_along(idx)) {
        pars <- as.numeric(params[i, ])
        qmin <- qgl(0, pars)
        if (is.finite(qmin)) {
            checkIdentical(0, pgl(qmin, pars))
            qout <- if (qmin < 0) 1.1 * qmin else 0.9 * qmin
            checkIdentical(0, pgl(qout, pars))
        }
    }

    ###################################
    # test values that are equal or larger than 1 quantile
    for (i in seq_along(idx)) {
        pars <- as.numeric(params[i, ])
        qmax <- qgl(1, pars)
        if (is.finite(qmax)) {
            checkIdentical(1, pgl(qmax, pars))
            qout <- if (qmax > 0) 1.1 * qmax else 0.9 * qmax
            checkIdentical(1, pgl(qout, pars))
        }
    }

    ###################################
    # test invalid options that should produce NaNs and warnings

    oo <- options(warn = -1)
    on.exit(options(oo))

    checkTrue(all(is.nan(pgl(c(0, 0, 0), iqr = 0))))
    checkTrue(all(is.nan(pgl(c(0, 0, 0), iqr = -1))))
    checkTrue(all(is.nan(pgl(c(0, 0, 0), chi = -1))))
    checkTrue(all(is.nan(pgl(c(0, 0, 0), chi = 1))))
    checkTrue(all(is.nan(pgl(c(0, 0, 0), chi = 2))))
    checkTrue(all(is.nan(pgl(c(0, 0, 0), chi = -2))))
    checkTrue(all(is.nan(pgl(c(0, 0, 0), xi = 0))))
    checkTrue(all(is.nan(pgl(c(0, 0, 0), xi = 1))))
    checkTrue(all(is.nan(pgl(c(0, 0, 0), xi = -1))))
    checkTrue(all(is.nan(pgl(c(0, 0, 0), xi = 2))))

    options(warn = 2)

    checkException(pgl(c(0, 0, 0), iqr = 0), silent = TRUE)
    checkException(pgl(c(0, 0, 0), iqr = -1), silent = TRUE)
    checkException(pgl(c(0, 0, 0), chi = -1), silent = TRUE)
    checkException(pgl(c(0, 0, 0), chi = 1), silent = TRUE)
    checkException(pgl(c(0, 0, 0), chi = 2), silent = TRUE)
    checkException(pgl(c(0, 0, 0), chi = -2), silent = TRUE)
    checkException(pgl(c(0, 0, 0), xi = 0), silent = TRUE)
    checkException(pgl(c(0, 0, 0), xi = 1), silent = TRUE)
    checkException(pgl(c(0, 0, 0), xi = -1), silent = TRUE)
    checkException(pgl(c(0, 0, 0), xi = 2), silent = TRUE)

    options(oo)

}
