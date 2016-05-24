
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# A copy of the GNU General Public License is available via WWW at
# http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
# writing to the Free Software Foundation, Inc., 59 Temple Place,
# Suite 330, Boston, MA  02111-1307  USA.

# Copyrights (C)
# for this R-port:
#   1999 - 2007, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:                 DESCRIPION:
#  tsTest                    Time Series Test Suite
# FUNCTION:                 DEPENDENCY TEST:
#  bdsTest                   Brock-Dechert-Scheinkman test for iid series
# FUNCTION:                 NONLINEARITY TESTS:
#  wnnTest                   White Neural Network Test for Nonlinearity
#  tnnTest                   Teraesvirta Neural Network Test for Nonlinearity
# FUNCTION:                 MORE TESTS ...
#  runsTest                  Runs test for detecting non-randomness [tseries]
################################################################################


################################################################################
# BUILTIN - PACKAGE DESCRIPTION:
#  Package: tseries
#  Version: 0.9-21
#  Date: 2004-04-23
#  Title: Time series analysis and computational finance
#  Author: Compiled by Adrian Trapletti <a.trapletti@bluewin.ch>
#  Maintainer: Kurt Hornik <Kurt.Hornik@R-project.org>
#  Description: Package for time series analysis and computational finance
#  Depends: R (>= 1.9.0), quadprog
#  License: GPL (see file COPYING)
#  Packaged: Thu Apr 22 16:32:16 2004; hornik
#  Notes: The runs.test is available as dependency test in the fBasics
#    Package runs.test = function (x)
#    Most of the functions are BUILTIN from Adrian Trapletti's R package
#    tseries
################################################################################


################################################################################
# REQUIREMENTS:             DESCRIPTION:
#  embed                     Required from fBasics.A0-SPlusCompatibility
################################################################################


tsTest =
function(x,
method = c("bds", "tnn", "wnn"), ...)
{   # A function implemented by Diethelm Wuertz

    # Load Library:
    # require(tseries)

    # Check Type:
    if (class(x) == "timeSeries") {
        if (dim(x)[2] > 1) stop("x must be an univariate time series")
    }
    x = as.vector(x)

    # Settings:
    method = match.arg(method)
    test = paste(method, "Test", sep = "")
    fun = match.fun(test)

    # Test:
    ans = fun(x = x, ...)

    # Return Value:
    ans

}


# ------------------------------------------------------------------------------


bdsTest =
function(x, m = 3, eps = NULL, title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Brock-Dechert-Scheinkman test for iid series

    # Notes:
    #   This function is a slightly modified copy of Adrian Trapletti's
    #   contributed function from his 'tseries' package.

    # FUNCTION:

    # Call:
    call = match.call()

    # Test
    test = list()

    # Data Set Name:
    DNAME = deparse(substitute(x))
    test$data.name = DNAME

    # Check Type:
    if (class(x) == "timeSeries") {
        if (dim(x)[2] > 1) stop("x must be an univariate time series")
    }
    x = as.vector(x)

    # Test:
    if (is.null(eps)) eps = seq(0.5*sd(x), 2*sd(x), length = 4)
    if (m < 2) stop("m is less than 2")
    if (length(eps) == 0) stop("invalid eps")
    if (any(eps <= 0)) stop("invalid eps")
    n = length(x)
    k = length(eps)
    cc = double(m+1)
    cstan = double(m+1)

    # Statistic:
    STATISTIC = NULL
    NAMES = NULL
    for(i in (1:k)) {
        res = .C("bdstest_main", as.integer(n), as.integer(m),
            as.double(x), as.double(cc), cstan = as.double(cstan),
            as.double(eps[i]), as.integer(0), PACKAGE = "fNonlinear")
        ans = res$cstan[2:m+1]
        STATISTIC = c(STATISTIC, ans)
        names.1 = rep(paste("eps[", i, "]", sep = ""), times = length(ans))
        names.2 = paste("m=", as.character(2:m), sep = "")
        NAMES = c(NAMES, paste(names.1, names.2))
    }
    # colnames(STATISTIC) = as.character(eps)
    # rownames(STATISTIC) = as.character(2:m)
    names(STATISTIC) = NAMES
    test$statistic = STATISTIC

    # P Value:
    PVAL = 2 * pnorm(-abs(STATISTIC))
    names(PVAL) = names(STATISTIC)
    # colnames(PVAL) = as.character(eps)
    # rownames(PVAL) = as.character(2:m)
    test$p.value = PVAL

    # METHOD = "BDS Test"

    PARAMETER = c(m, eps)
    names(PARAMETER) = c(
        "Max Embedding Dimension",
        paste("eps[", 1:length(eps), "]", sep = "") )
    test$parameter = PARAMETER

    # Add:
    if (is.null(title)) title = "BDS Test"
    if (is.null(description)) description = description()

    # Return Value:
    new("fHTEST",
        call = call,
        data = list(x = x),
        test = test,
        title = as.character(title),
        description = as.character(description) )
}


# ------------------------------------------------------------------------------


wnnTest =
function(x, lag = 1, qstar = 2, q = 10, range = 4,
title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   White's Neural Network Test for Nonlinearity

    # Notes:
    #   This function is a slightly modified copy of Adrian Trapletti's
    #   contributed function from his 'tseries' package.

    # FUNCTION:

    # Call:
    CALL = match.call()

    # Test
    test = list()

    # Data Set Name:
    DNAME = deparse(substitute(x))
    test$data.name = DNAME

    # Check Type:
    if (class(x) == "timeSeries") {
        if (dim(x)[2] > 1) stop("x must be an univariate time series")
    }
    x = as.vector(x)

    # Test - This part comes from A. Trapletti's Code:
    if (lag < 1) stop("minimum lag is 1")
    t = length(x)
    y = embed(x, lag + 1)
    xnam = paste("y[,", 2:(lag+1), "]", sep = "")
    fmla = as.formula(paste("y[,1]~", paste(xnam, collapse = "+")))
    rr = lm(fmla)
    u = residuals(rr)
    ssr0 = sum(u^2)
    max = range/2
    gamma = matrix(runif ((lag+1)*q, -max, max), lag+1, q)
    phantom = (1 + exp(-(cbind(rep(1, t-lag), y[, 2:(lag+1)]) %*% gamma)))^(-1)
    # Changed to be compatible with SPlus:
    # phantomstar = as.matrix(prcomp(phantom, scale = TRUE)$x[, 2:(qstar+1)])
    phantomstar = as.matrix(prcomp(scale(phantom))$x[, 2:(qstar+1)])
    xnam2 = paste("phantomstar[,", 1:qstar, "]", sep = "")
    xnam2 = paste(xnam2, collapse = "+")
    fmla = as.formula(paste("u~", paste(paste(xnam, collapse = "+"),
        xnam2, sep = "+")))
    rr = lm(fmla)
    v = residuals(rr)
    ssr = sum(v^2)

    # Statistic:
    STATISTIC1 = t * log(ssr0/ssr)
    STATISTIC2 = ((ssr0-ssr)/qstar)/(ssr/(t-lag-qstar))
    STATISTIC = c(STATISTIC1, STATISTIC2)
    names(STATISTIC) = c("Chi-squared", "F")
    test$statistic = STATISTIC

    # P Values:
    PVAL1 = 1 - pchisq(STATISTIC1, qstar)
    PVAL2 = 1 - pf(STATISTIC2, qstar, t-lag-qstar)
    PVAL = c(PVAL1, PVAL2)
    names(PVAL) = c("Chi-squared", "F")
    test$p.value = PVAL

    # Parameter:
    PARAMETER = c(lag, q, range, qstar, t-lag-qstar)
    names(PARAMETER) = c("lag", "q", "range", "qstar|df", "t-lag-qstar|df")
    test$parameter = PARAMETER

    # Add:
    if (is.null(title)) title = "White Neural Network Test"
    if (is.null(description)) description = description()

    # Return Value:
    new("fHTEST",
        call = CALL,
        data = list(x = x, y = y),
        test = test,
        title = as.character(title),
        description = as.character(description) )
}


# ------------------------------------------------------------------------------


tnnTest =
function(x, lag = 1, title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Teraesvirta's Neural Network Test for Nonlinearity

    # Notes:
    #   This function is a slightly modified copy of Adrian Trapletti's
    #   contributed function from his 'tseries' package.

    # FUNCTION:

    # Call:
    call = match.call()

    # Test
    test = list()

    # Data Set Name:
    DNAME = deparse(substitute(x))
    test$data.name = DNAME

    # Check Type:
    if (class(x) == "timeSeries") {
        if (dim(x)[2] > 1) stop("x must be an univariate time series")
    }
    x = as.vector(x)

    # Test - This part comes from A. Trapletti's Code:
    if (lag < 1) stop("minimum lag is 1")
    t = length(x)
    y = embed(x, lag+1)
    xnam = paste("y[,", 2:(lag+1), "]", sep = "")
    fmla = as.formula(paste("y[,1]~", paste(xnam, collapse = "+")))
    rr = lm(fmla)
    u = residuals(rr)
    ssr0 = sum(u^2)
    xnam2 = NULL
    m = 0
    for (i in (1:lag)) {
        for (j in (i:lag)) {
            xnam2 = c(xnam2, paste("I(y[,",i+1, "]*y[, ",j+1,"])", sep = ""))
            m = m+1
        }
    }
    xnam2 = paste(xnam2, collapse="+")
    xnam3 = NULL
    for (i in (1:lag)) {
        for (j in (i:lag)) {
            for (k in (j:lag)) {
                xnam3 = c(xnam3,
                paste("I(y[,", i+1, "]*y[,", j+1, "]*y[,", k+1, "])", sep = ""))
                m = m+1
            }
        }
    }
    xnam3 = paste(xnam3,collapse="+")
    fmla = as.formula(paste("u~", paste(paste(xnam, collapse = "+"),
        xnam2, xnam3, sep = "+")))
    rr = lm(fmla)
    v = residuals(rr)
    ssr = sum(v^2)

    #Statistic:
    STATISTIC1 = t*log(ssr0/ssr)
    STATISTIC2 = ((ssr0-ssr)/m)/(ssr/(t-lag-m))
    STATISTIC = c(STATISTIC1, STATISTIC2)
    names(STATISTIC) = c("Chi-squared", "F")
    test$statistic = STATISTIC

    # P Value:
    PVAL1 = 1 - pchisq(STATISTIC1, m)
    PVAL2 = 1 - pf(STATISTIC2, m, t-lag-m)
    PVAL = c(PVAL1, PVAL2)
    names(PVAL) = c("Chi-squared", "F")
    test$p.value = PVAL

    # PARAMETER:
    PARAMETER = c(lag, m, t-lag-m)
    names(PARAMETER) = c("lag", "m|df", "t-lag-m|df")
    test$parameter = PARAMETER

    # Add:
    if (is.null(title)) title = "Teraesvirta Neural Network Test"
    if (is.null(description)) description = description()

    # Return Value:
    new("fHTEST",
        call = call,
        data = list(x = x),
        test = test,
        title = as.character(title),
        description = as.character(description) )
}


################################################################################


runsTest =
function(x)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Performs a runs test

    # Arguments:
    #   x - a numeric vector of data values.

    # Notes:
    #   Implementing Trapletti's tseries R-Package

    # Note:
    #   We consider the signs of x in the series, the zeros will be
    #   discarded. In addition we have to factor the data for runs.test().

    # FUNCTION:

    # Convert Type:
    if (class(x) == "fREG") x = residuals(x)
    x = as.vector(x)

    # runs.test() copied from A. Traplettis tseries package
    runs.test =
    function (x, alternative = c("two.sided", "less", "greater"))
    {
        if (!is.factor(x)) stop("x is not a factor")
        if (any(is.na(x))) stop("NAs in x")
        if (length(levels(x)) != 2) stop("x does not contain dichotomous data")
        alternative = match.arg(alternative)
        DNAME = deparse(substitute(x))
        n = length(x)
        R = 1 + sum(as.numeric(x[-1] != x[-n]))
        n1 = sum(levels(x)[1] == x)
        n2 = sum(levels(x)[2] == x)
        m = 1 + 2 * n1 * n2/(n1 + n2)
        s = sqrt(2 * n1 * n2 * (2 * n1 * n2 - n1 - n2)/((n1 + n2)^2 *
            (n1 + n2 - 1)))
        STATISTIC = (R - m)/s
        METHOD = "Runs Test"
        if (alternative == "two.sided")
            PVAL = 2 * pnorm(-abs(STATISTIC))
        else if (alternative == "less")
            PVAL = pnorm(STATISTIC)
        else if (alternative == "greater")
            PVAL = pnorm(STATISTIC, lower.tail = FALSE)
        else stop("irregular alternative")
        names(STATISTIC) = "Standard Normal"
        structure(list(
            statistic = STATISTIC,
            alternative = alternative,
            p.value = PVAL,
            method = METHOD,
            data.name = DNAME),
            class = "htest") }

    # Result:
    x = sign(x)
    x = x[x != 0]
    x = factor(x)
    ans = runs.test(x = x)

    # Return Value:
    ans
}


################################################################################

