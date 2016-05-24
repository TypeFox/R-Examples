
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
# FUNCTION:                ADF TESTS:
#  adfTest                  ADF unit root test using Banarjee's test statistics
#  unitrootTest             ADF unit root test using McKinnon's test statistics
# FUNCTION:                UNITROOT TEST SUITE:
#  .urTest                  Unit Root Test Suite
################################################################################



adfTest =
function(x, lags = 1, type = c("nc", "c", "ct"), title = NULL,
description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Tests the null hypothesis of a unit root in y.

    # Arguments:
    #   x - numeric vector
    #   type - specifies the regression model to be estimatied and the
    #       null hypothesis, "nc" no constant and no trend, "c" add
    #       constant, "ct" add constant and trend.
    #   lags - specifies the number of lagged differences of x to be
    #       included in the regression model. If 'lags' = h, a term
    #       sum_{i=1}^{h-1} beta_i * diff(x)_(t-h) is added to the
    #       regression equation.

    # Value:
    #   A list of class "htest" containing the following components:
    #   statistic - the value of the test statistic (t-statistic)
    #   parameter - the number of lags.
    #   p.value - the p-value of the test
    #   method - a character string indicating what type of test was performed
    #   data.name - a character string giving the name of the data y

    # Reference:
    #   S. E. SAID and D. A. DICKEY (1984): Testing for Unit Roots in
    #   Autoregressive-Moving Average Models of Unlag.diffnown Order.
    #   Biometrika 71, 599–607.

    # Source:
    #   This function is an augmented version of Adrian Trapletti's
    #   function adf.test() which considers type "ct" only. We have added
    #   the types "c" and "nc" together with the appropriate statistics.

    # Call:
    CALL = match.call()

    # Test:
    test = list()

    # Data Set Name:
    DNAME = deparse(substitute(x))
    test$data.name = DNAME

    # Transform:
    if (class(x) == "timeSeries") x = series(x)
    x = as.vector(x)

    # Check Arguments:
    if (lags < 0) stop("Lags are negative")

    # Settings:
    doprint = FALSE
    type = type[1]
    lags = lags + 1
    y = diff(x)
    n = length(y)
    z = embed(y, lags)
    y.diff = z[, 1]
    y.lag.1 = x[lags:n]
    tt = lags:n

    # Regression:
    if (lags > 1) {
        y.diff.lag = z[,2:lags]
        if (type == "nc"){
            res = lm(y.diff ~ y.lag.1 - 1 + y.diff.lag) }
        if (type == "c"){
            res = lm(y.diff ~ y.lag.1 + 1 +  y.diff.lag) }
        if (type == "ct") {
            res = lm(y.diff ~ y.lag.1 + 1 + tt + y.diff.lag) }
    } else {
        if (type == "nc") {
            res = lm(y.diff ~ y.lag.1 - 1) }
        if (type == "c"){
            res = lm(y.diff ~ y.lag.1 + 1) }
        if (type == "ct") {
            res = lm(y.diff ~ y.lag.1 + 1  + tt)  }
    }

    # Regression Summary:
    res.sum = summary(res)
    if (doprint) print(res.sum)

    # Statistic:
    if (type == "nc") coefNum = 1 else coefNum = 2
    STAT = res.sum$coefficients[coefNum, 1] / res.sum$coefficients[coefNum, 2]
    names(STAT) = "Dickey-Fuller"
    test$statistic = STAT

    # P Value:
    if (type == "nc")
        table = cbind(
            c(-2.66, -2.26, -1.95, -1.60, +0.92, +1.33, +1.70, +2.16),
            c(-2.62, -2.25, -1.95, -1.61, +0.91, +1.31, +1.66, +2.08),
            c(-2.60, -2.24, -1.95, -1.61, +0.90, +1.29, +1.64, +2.03),
            c(-2.58, -2.23, -1.95, -1.62, +0.89, +1.29, +1.63, +2.01),
            c(-2.58, -2.23, -1.95, -1.62, +0.89, +1.28, +1.62, +2.00),
            c(-2.58, -2.23, -1.95, -1.62, +0.89, +1.28, +1.62, +2.00))
    if (type == "c")
        table = cbind(
            c(-3.75, -3.33, -3.00, -2.63, -0.37, +0.00, +0.34, +0.72),
            c(-3.58, -3.22, -2.93, -2.60, -0.40, -0.03, +0.29, +0.66),
            c(-3.51, -3.17, -2.89, -2.58, -0.42, -0.05, +0.26, +0.63),
            c(-3.46, -3.14, -2.88, -2.57, -0.42, -0.06, +0.24, +0.62),
            c(-3.44, -3.13, -2.87, -2.57, -0.43, -0.07, +0.24, +0.61),
            c(-3.43, -3.12, -2.86, -2.57, -0.44, -0.07, +0.23, +0.60))
    if (type == "ct")
        table = cbind(
            c(-4.38, -3.95, -3.60, -3.24, -1.14, -0.80, -0.50, -0.15),
            c(-4.15, -3.80, -3.50, -3.18, -1.19, -0.87, -0.58, -0.24),
            c(-4.04, -3.73, -3.45, -3.15, -1.22, -0.90, -0.62, -0.28),
            c(-3.99, -3.69, -3.43, -3.13, -1.23, -0.92, -0.64, -0.31),
            c(-3.98, -3.68, -3.42, -3.13, -1.24, -0.93, -0.65, -0.32),
            c(-3.96, -3.66, -3.41, -3.12, -1.25, -0.94, -0.66, -0.33))
    table = t(table)
    tablen = dim(table)[2]
    tableT = c(25, 50, 100, 250, 500, 1e+05)
    tablep = c(0.01, 0.025, 0.05, 0.1, 0.9, 0.95, 0.975, 0.99)
    tableipl = numeric(tablen)
    for (i in (1:tablen))
        tableipl[i] = approx(tableT, table[, i], n, rule = 2)$y
    PVAL = approx(tableipl, tablep, STAT, rule = 2)$y
    if (is.na(approx(tableipl, tablep, STAT, rule = 1)$y)) {
        if (PVAL == min(tablep)) {
            warning("p-value smaller than printed p-value")
        } else {
            warning("p-value greater than printed p-value")
        }
    }
    names(PVAL) = ""
    test$p.value = PVAL

    # Parameter:
    PARAMETER = lags - 1
    names(PARAMETER) = "Lag Order"
    test$parameter = PARAMETER

    # Add:
    if (is.null(title)) title = "Augmented Dickey-Fuller Test"
    if (is.null(description)) description = date()

    # Add Regression:
    test$lm = res

    # Return Value:
    new("fHTEST",
        call = CALL,
        data = list(x = x),
        test = test,
        title = as.character(title),
        description = description()
        )
}


# ------------------------------------------------------------------------------


unitrootTest =
function(x, lags = 1, type = c("nc", "c", "ct"), title = NULL,
description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Tests the null hypothesis of a unit root in x.

    # Arguments:
    #   x - numeric vector
    #   type - specifies the regression model to be estimatied and the
    #       null hypothesis, "nc" no constant and no trend, "c" add
    #       constant, "ct" add constant and trend.
    #   lags - specifies the number of lagged differences of x to be
    #       included in the regression model. If 'lags' = h, a term
    #       sum_{i=1}^{h-1} beta_i * diff(x)_(t-h) is added to the
    #       regression equation.

    # Value:
    #   A list with class "htest" containing the following components:
    #   statistic - the value of the test statistic (t-statistic)
    #   parameter - the number of lags.
    #   p.value - the p-value of the test
    #   method - a character string indicating what "trend" type of
    #       the test was performed
    #   data.name - a character string giving the name of the data y

    # Reference:
    #   Said S.E., Dickey D.A. (1984): Testing for Unit Roots in
    #   Autoregressive-Moving Average Models of Unlag.diffnown Order.
    #   Biometrika 71, 599–-607.

    # Source:
    #   This function is an augmented version of Adrian Trapletti's
    #   function adf.test() which considers trend "ct" only. We have added
    #   the trend types "c" and "nc" together with the appropriate statistics.

    # FUNCTION:

    # Call:
    CALL = match.call()

    # Test:
    test = list()

    # Data Set Name:
    DNAME = deparse(substitute(x))
    test$data.name = DNAME

    # Transform:
    if (class(x) == "timeSeries") x = series(x)
    x = as.vector(x)

    # Check Arguments:
    if (lags < 0) stop("Lags are negative")

    # Settings:
    type = type[1]
    lags = lags + 1
    y = diff(x)
    n = length(y)
    z = embed(y, lags)
    y.diff = z[, 1]
    y.lag.1 = x[lags:n]
    tt = lags:n

    # Regression:
    if (lags > 1) {
        y.diff.lag = z[,2:lags]
        if (type == "nc"){
            res = lm(y.diff ~ y.lag.1 - 1 + y.diff.lag) }
        if (type == "c"){
            res = lm(y.diff ~ y.lag.1 + 1 +  y.diff.lag) }
        if (type == "ct") {
            res = lm(y.diff ~ y.lag.1 + 1 + tt + y.diff.lag) }
        if (type == "ctt") {
            res = lm(y.diff ~ y.lag.1 + 1 + tt + tt^2 + y.diff.lag) }
    } else {
        if (type == "nc") {
            res = lm(y.diff ~ y.lag.1 - 1) }
        if (type == "c"){
            res = lm(y.diff ~  y.lag.1 + 1) }
        if (type == "ct") {
            res = lm(y.diff ~ y.lag.1 + 1  + tt)  }
        if (type == "ctt") {
            res = lm(y.diff ~ y.lag.1 + 1 + tt + tt^2) }
    }
    res.sum = summary(res)
    test$regression = res.sum

    # Statistic:
    if (type == "nc") coefNum = 1 else coefNum = 2
    STATISTIC =
        res.sum$coefficients[coefNum, 1] / res.sum$coefficients[coefNum, 2]
    names(STATISTIC) = "DF"
    test$statistic = STATISTIC

    # P Value:
    if (type == "nc") { itv = 1 }
    if (type == "c")  { itv = 2 }
    if (type == "ct") { itv = 3 }
    if (type == "ctt"){ itv = 4 }
    # Statistic == "t" : itt = 1
    PVAL1 =
        .urcval(arg = STATISTIC, nobs = n, niv = 1, itt = 1, itv = itv, nc = 2)
    # Statistic == "n" : itt = 2
    PVAL2 =
        .urcval(arg = STATISTIC, nobs = n, niv = 1, itt = 2, itv = itv, nc = 2)
    PVAL = c(PVAL1, PVAL2)
    names(PVAL) = c("t", "n")
    test$p.value = PVAL

    # Parameter:
    PARAMETER = lags - 1
    names(PARAMETER) = "Lag Order"
    test$parameter = PARAMETER

    # Add:
    if (is.null(title)) title = "Augmented Dickey-Fuller Test"
    if (is.null(description)) description = date()

    # Return Value:
    new("fHTEST",
        call = CALL,
        data = list(x = x),
        test = test,
        title = as.character(title),
        description = description()
        )
}


################################################################################


.urTest =
function(x, method = c("unitroot", "adf", "urers", "urkpss", "urpp",
"ursp", "urza"), title = NULL, description = NULL, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Unit Root Test Suite

    # FUNCTION:

    # Match Function:
    funTest = match.fun(paste(method[1], "Test", sep = ""))

    # Test:
    ans = funTest(x = x, ...)

    # Add:
    if (!is.null(title)) ans@title = as.character(title)
    if (!is.null(description)) ans@description = description()

    # Return Value:
    ans
}


################################################################################

