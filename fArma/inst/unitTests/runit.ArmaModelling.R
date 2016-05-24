
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA

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
# FUNCTION:               SIMULATION AND FITTING:
#  'fARMA'                 S4 Class representation for "fARMA" objects
#  armaSim                 Simulates an ARIMA time series process
#  armaFit                 Fits parameters for ARMA Time Series process
# S3 METHOD:              PREDICTION:
#  predict.fARMA           S3: Predicts from an ARMA time series prrocess
# GENERIC METHODS:        PRINT - PLOT - SUMMARY METHODS:
#  show.fARMA              S4: Prints a fitted ARMA time series object
#  plot.fARMA              S3: Plots stylized facts of a fitted ARMA object
#  summary.fARMA           S3: Summarizes a fitted ARMA time series object
# S3 METHOD:              ADDON METHODS:
#  coef.fARMA              S3: Returns coefficidents from a fitted ARMA object
#  coefficients.fARMA      S3: Synonyme for coef.fARMA
#  fitted.fARMA            S3: Returns fitted values from a fitted ARMA object
#  residuals.fARMA         S3: Returns residuals from a fitted ARMA object
################################################################################


test.fARMA =
function()
{
    # Slot Names:
    Names = getSlots("fARMA")
    target = names(Names)
    current = c(
        "call",
        "formula",
        "method",
        "parameter",
        "data",
        "fit",
        "residuals",
        "fitted",
        "title",
        "description")
    checkIdentical(target, current)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.armaSim =
function()
{
    # armaSim(model = list(ar = c(0.5, -0.5), d = 0, ma = 0.1), n = 100,
    #   positions = NULL, innov = NULL, n.start = 100, start.innov = NULL,
    #   rand.gen = rnorm, rseed = NULL, addControl = FALSE, ...)

    # ts: Simulate ARMA(2,1):
    # Note, if "positions=NULL",
    #    then a 'ts' object will be returned ...
    ts = armaSim(n = 25)
    class(ts)
    print(ts)
    ts = armaSim(n = 25, addControl = TRUE)
    print(ts)

    # timeSeries: Simulate ARMA(2,1):
    # Note, if "positions" is a timeDate object,
    #   then a 'timeSeries' object will be returned ...
    tS = armaSim(n = 12)
    class(tS)
    print(tS)
    tS = armaSim(n = 12, addControl = TRUE)
    print(tS)

    # ts: Simulate t4-ARMA(2,1):
    ts = armaSim(n = 25, rand.gen = rt, df = 4, rseed = 4711)
    print(ts)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.ar2Fit =
function()
{
    # Simulate AR(2):
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = armaSim(model = list(ar = c(0.5, -0.5)), n = 1000)
    head(x)

    # method = c("mle", "ols")
    args(ar)

    # AR(2) - OLS Fit:
    # R's internal function ar() will be used ...
    object = armaFit(formula = ~ ar(2), data = x, method = "ols")
    print(object)
    target = as.vector(round(coef(object), 1))
    current = c(0.5, -0.5, 0)
    checkEqualsNumeric(target, current)

    # AR(2) - MLE Fit:
    # R's internal function ar() will be used ...
    object = armaFit(formula = ~ ar(2), data = x, method = "mle")
    print(object)
    target = as.vector(round(coef(object), 1))
    current = c(0.5, -0.5, 0)
    checkEqualsNumeric(target, current)

    # For the expert ...
    # Note, also other methods can be used supported by ar():

    # Yule-Walker Fit:
    object = armaFit(formula = ~ ar(2), data = x, method = "yw")
    print(object)
    target = as.vector(round(coef(object), 1))
    current = c(0.5, -0.5, 0)
    checkEqualsNumeric(target, current)

    # Burg 1 Fit:
    object = armaFit(formula = ~ ar(2), data = x, method = "burg1")
    print(object)
    target = as.vector(round(coef(object), 1))
    current = c(0.5, -0.5, 0)
    checkEqualsNumeric(target, current)

    # Burg 2 Fit:
    object = armaFit(formula = x ~ ar(2), data = x, method = "burg2")
    print(object)
    target = as.vector(round(coef(object), 1))
    current = c(0.5, -0.5, 0)
    checkEqualsNumeric(target, current)

    # Note, also arma() or arima() formulas can be applied:
    # In these cases R's internal function arima() will be used ...
    args(arima)

    # CSS-ML Fit:
    object = armaFit(formula = ~ arima(2, 0, 0), data =  x, method = "mle")
    print(object)
    target = as.vector(round(coef(object), 1))
    current = c(0.5, -0.5, 0)
    checkEqualsNumeric(target, current)

    # Alternatively you can use also method=CSS-ML Fit:
    object = armaFit(formula = ~ arima(2, 0, 0), data =  x, method = "CSS-ML")
    print(object)
    target = as.vector(round(coef(object), 1))
    current = c(0.5, -0.5, 0)
    checkEqualsNumeric(target, current)

    # or method=CSS Fit:
    object = armaFit(formula = ~ arima(2, 0, 0), data = x, method = "CSS")
    print(object)
    target = as.vector(round(coef(object), 1))
    current = c(0.5, -0.5, 0)
    checkEqualsNumeric(target, current)

    # or method=ML Fit:
    object = armaFit(formula = ~ arima(2, 0, 0), data = x, method = "ML")
    print(object)
    target = as.vector(round(coef(object), 1))
    current = c(0.5, -0.5, 0)
    checkEqualsNumeric(target, current)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.ar2Report =
function()
{
    # Simulate AR(2):
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = armaSim(model = list(ar = c(0.5, -0.5)), n = 1000)
    head(x)

    # Fit:
    object = armaFit(formula = ~ ar(2), data =  x, method = "mle")

    # Short Print Report:
    print(object)

    # Plot: Standardized Residuals, ACF, QQ-Plot, Ljung-Box p-Values
    par(mfrow = c(2, 2), cex = 0.7)
    plot(object, which = "all")

    # Try:
    par(mfrow = c(1, 1))
    plot(object, which = 1)
    plot(object, which = 2)
    plot(object, which = 3)
    plot(object, which = 4)

    # Interactive Plot:
    # plot(object)

    # Summary Method - No Plot:
    summary(object, doplot = FALSE)

    # Summary Method - Including All Plots:
    par(mfrow = c(2, 2), cex = 0.7)
    summary(object, doplot = TRUE, which = "all")

    # Extractor Functions - Get Values:
    coefficients(object)
    coef(object)
    fitted = fitted(object)
    class(fitted)
    head(fitted)
    residuals = residuals(object)
    class(residuals)
    head(residuals)

    # Should we have ?
    # getCoef
    # getResiduals
    # getFitted

    # Predict Method:
    # predict.fARMA is now in namespace as an hidden method
    # args(predict.fARMA)
    predict(object)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.ma2Fit =
function()
{
    # Simulate MA(2):
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = armaSim(model = list(d = 0, ma = c(0.5, -0.5)), n = 5000)

    # To Fit a MA Model use ma(q), arma(0,q) or arima(0, 0, q):
    object = armaFit(formula = ~ ma(2), data = x)
    print(object)
    target = as.vector(round(coef(object), 1))
    current = c(0.5, -0.5, 0)
    checkEqualsNumeric(target, current)

    # Note, also arma() or arima() formulas can be applied:

    # CSS-ML Fit:
    object = armaFit(formula = ~ arima(0, 0, 2), data = x, method = "CSS-ML")
    print(object)
    target = as.vector(round(coef(object), 1))
    current = c(0.5, -0.5, 0)
    checkEqualsNumeric(target, current)

    # CSS Fit:
    object = armaFit(formula = ~ arima(0, 0, 2), data = x, method = "CSS")
    print(object)
    target = as.vector(round(coef(object), 1))
    current = c(0.5, -0.5, 0)
    checkEqualsNumeric(target, current)

    # ML fit:
    object = armaFit(formula = ~ arima(0, 0, 2), data = x, method = "ML")
    print(object)
    target = as.vector(round(coef(object), 1))
    current = c(0.5, -0.5, 0)
    checkEqualsNumeric(target, current)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.ma2Report =
function()
{
    # Simulate MA(2):
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = armaSim(model = list(d = 0, ma = c(0.5, -0.5)), n = 5000)

    # To Fit a MA Model use ma(q), arma(0,q) or arima(0, 0, q):
    object = armaFit(formula = ~ ma(2), data = x)

    # Report:
    print(object)

    # Plot: Standardized Residuals, ACF, QQ-Plot, Ljung-Box p-Values
    par(mfrow = c(2, 2), cex = 0.7)
    plot(object, which = "all")

    # Summary:
    summary(object, doplot = FALSE)

    # Get Values:
    coefficients(object)
    coef(object)
    fitted(object)[1:10]
    residuals(object)[1:10]

    # Predict:
    predict(object)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.arma21Fit =
function()
{
    # Simulate ARMA(2,1):
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = armaSim(model = list(ar = c(0.5, -0.5), ma = 0.1), n = 1000)

    # Fit:
    object = armaFit(formula = ~ arma(2, 1), data =  x, method = "mle")
    print(object)
    target = as.vector(round(coef(object), 1))
    print(target)
    current = c(0.5, -0.5, 0.1, 0)
    checkEqualsNumeric(target, current)

    # Note, also arima() formulas can be applied:

    # "mle" == "ML" Fit:
    object = armaFit(formula = ~ arima(2, 0, 1), data =  x)
    print(object)
    target = as.vector(round(coef(object), 1))
    print(target)
    current = c(0.5, -0.5, 0.1, 0)
    checkEqualsNumeric(target, current)

    # CSS Fit:
    object = armaFit(formula = ~ arima(2, 0, 1), data =  x, method = "CSS")
    print(object)
    target = as.vector(round(coef(object), 1))
    print(target)
    current = c(0.5, -0.5, 0.1, 0)
    checkEqualsNumeric(target, current)

    # CSS-ML Fit:
    object = armaFit(formula = ~ arima(2, 0, 1), data =  x, method = "CSS-ML")
    print(object)
    target = as.vector(round(coef(object), 1))
    print(target)
    current = c(0.5, -0.5, 0.1, 0)
    checkEqualsNumeric(target, current)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.arma21Report =
function()
{
    # Simulate ARMA(2, 1):
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = armaSim(model = list(ar = c(0.5, -0.5), ma = 0.1), n = 1000)

    # Fit:
    object = armaFit(formula = ~ arma(2, 1), data = x)

    # Report:
    print(object)

    # Plot:
    par(mfrow = c(2, 2), cex = 0.7)
    plot(object, which = "all")

    # Summary:
    summary(object, doplot = FALSE)

    # Get Values:
    coefficients(object)
    coef(object)
    fitted(object)[1:10]
    residuals(object)[1:10]

    # Predict:
    predict(object)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.arima211Fit =
function()
{
    # Simulate ARIMA(2, 1, 1):
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = armaSim(model = list(ar = c(0.5, -0.5), d = 1, ma = 0.1), n = 1000)

    # CSS-ML Fit:
    object = armaFit(formula = ~ arima(2, 1, 1), data = x, method = "CSS-ML")
    print(object)
    target = as.vector(round(coef(object), 1))
    print(target)
    current = c(0.5, -0.5, 0.1)
    checkEqualsNumeric(target, current)

    # CSS Fit:
    object = armaFit(formula = ~ arima(2, 1, 1), data = x, method = "CSS")
    print(object)
    target = as.vector(round(coef(object), 1))
    print(target)
    current = c(0.5, -0.5, 0.1)
    checkEqualsNumeric(target, current)

    # ML Fit:
    object = armaFit(formula = ~ arima(2, 1, 1), data = x, method = "ML")
    print(object)
    target = as.vector(round(coef(object), 1))
    print(target)
    current = c(0.5, -0.5, 0.1)
    checkEqualsNumeric(target, current)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.arima211Report =
function()
{
    # Simulate:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = armaSim(model = list(ar = c(0.5, -0.5), d = 1, ma = 0.1), n = 1000)

    # mle Integrated ARMA Fit:
    object = armaFit(formula = ~ arima(2, 1, 1), data = x)

    # Report:
    print(object)

    # Plot:
    par(mfrow = c(2, 2), cex = 0.7)
    plot(object, which = "all")

    # Summary
    summary(object, doplot = FALSE)

    # Get Values:
    coefficients(object)
    coef(object)
    fitted(object)[1:10]
    residuals(object)[1:10]

    # Predict:
    predict(object)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.arfima00Fit =
function()
{
    # Simulate ARFIMA(0, 0):
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = armaSim(model = list(d = 0.3), n = 1000)

    # Fit:
    object = armaFit(formula = ~ arfima(0, 0), data = x)
    print(object)
    target = as.vector(round(coef(object), 1))
    print(target)
    current = 0.3
    checkEqualsNumeric(target, current)

    # Parameter:
    target = unlist(object@parameter)
    print(target)
    current = c(include.mean = 1, M = 100, h = -1)
    checkIdentical(target, current)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.arfima00Report =
function()
{
    # Simulate:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = armaSim(model = list(d = 0.3), n = 1000)

    # Fit:
    object = armaFit(formula = ~ arfima(0, 0), data = x, M = 50, h = -1)

    # Report:
    print(object)

    # Plot:
    # plot(object, which = "all")         # CHECK not yet implemented

    # Summary:
    summary(object, doplot = FALSE)       # CHECK use always doplot=FALSE

    # Get Values:
    coefficients(object)
    coef(object)
    fitted = fitted(object)
    class(fitted)
    tail(fitted)                          # CHECK head yields NA's
    residuals = residuals(object)
    class(residuals)
    tail(residuals)

    # Predict:
    # predict(object)                     # CHECK not yet implemented

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.armaFormula =
function()
{
    # This example shows variants of alternative usage
    #   of the argument formula ...

    # Load From Ecofin Package:
    TS = MSFT
    head(TS)
    class(TS)
    colnames(TS)

    # Fit:
    armaFit(formula = diff(log(Close)) ~ ar(5), data = TS)

    # Fit:
    # Note, data may be a timeSeries object ...
    armaFit(Close ~ ar(5), data = returns(TS, digits = 12))

    # Fit:
    TS.RET = returns(TS, digits = 12)
    # Note, data may be a timeSeries object ...
    armaFit(Close ~ ar(5), TS.RET)

    # Fit:
    # Note, data may be a 'data.frame' ...
    armaFit(Close ~ ar(5), as.data.frame(TS.RET))

    # Fit:
    # Note, data may be a 'numeric' vector ...
    armaFit(~ ar(5), as.vector(TS.RET[, "Close"]))

    # Fit:
    # Note, data may be an object of class 'ts' ...
    armaFit(~ ar(5), as.ts(TS.RET)[, "Close"])

    # Fit:
    TS.RET = returns(TS)
    colnames(TS.RET)
    # attach(TS.RET)                    # CHECK doesn't work under RUnit ...
    attach(TS.RET)
    head(Close)
    armaFit(Close ~ ar(5))

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.armaArguments =
function()
{
    # armaFit(
    #   formula, data, method = c("mle", "ols"), include.mean = TRUE,
    #   fixed = NULL, title = NULL, description = NULL, ...)

    # arima(
    #   x, order = c(0, 0, 0), seasonal = list(order = c(0, 0, 0), period = NA),
    #   xreg = NULL, include.mean = TRUE, transform.pars = TRUE,
    #   fixed = NULL, init = NULL, method = c("CSS-ML", "ML", "CSS"),
    #   n.cond, optim.control = list(), kappa = 1e+06)

    # Simulate:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = armaSim(model = list(ar = c(0.5, -0.5), d = 0, ma = 0.1), n = 1000)

    # Include Mean - by default:
    armaFit(~ arma(2, 1), data = x)

    # Don't include Mean:
    armaFit(~ arma(2, 1), data = x, include.mean = FALSE)

    # Full Model:
    armaFit(~ arma(2, 1), data = x)

    # Fixed - AR(2[2]) Subset Model:
    #                                         ar1 ar2 ma1 intercept
    armaFit(~ arma(2, 1), data = x, fixed = c(0.5, NA, NA, NA))

    # Return Value:
    return()
}


################################################################################

