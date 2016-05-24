
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
# S3 METHOD:              PREDICTION:
#  predict.fARMA           S3: Predicts from an ARMA time series prrocess
#  .arPpredict              Internal function called by predict.fARMA
#  .arimaPpredict           Internal function called by predict.fARMA
#  .arfimaPredict           Internal function - Not yet implemented
################################################################################


predict.fARMA =
function (object, n.ahead = 10, n.back = 50, conf = c(80, 95),
doplot = TRUE, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Predicts from an ARMA Time Series Process

    # Example:
    #   x = armaSim(n = 500)
    #   object = armaFit(formula = x ~ arima(2, 0, 1))
    #   predict(object)

    # FUNCTION:

###     # OX Arfima:
###     if (object@call[[1]] == "arfimaOxFit") {
###         # .arfimaOxPredict(object, n.ahead = 10, n.back = 50, trace = FALSE)
###         ans = .arfimaOxPredict(object, n.ahead, n.back, ...)
###         return(ans)
###     }

    # Predict "ar":
    if (object@fit$tsmodel == "ar") {
        pred = .arPredict(object, n.ahead, se.fit = TRUE, ...)
    }

    # Predict "arima":
    if (object@fit$tsmodel == "arima") {
        pred = .arimaPredict(object, n.ahead, se.fit = TRUE, ...)
    }

    # Predict "arfima":
    if (object@fit$tsmodel == "arfima") {
        warning(" Prediction for ARFIMA not yet implemented")
        return()
    }

###     # Predict "arfima" from Ox:
###     if (object@fit$tsmodel == "arfimaOX") {
###         pred = .arfimaOxPredict(object, n.ahead, ...)
###     }

    # Prediction:
    names(pred$pred) = names(pred$se) = NULL
    ans = list(pred = pred$pred, se = pred$se)

    # Plot:
    if (doplot) {

        # Data:
        data = as.ts(object@data$x)
        freq = frequency(data)
        start = start(data)
        n = length(data)

        # Fit Slot:
        options(warn = -1)
        fit = object@fit
        class(fit) = fit$class

        # Upper and Lower Bands:
        nint = length(conf)
        upper = lower = matrix(NA, ncol = nint, nrow = length(pred$pred))
        for (i in 1:nint) {
            qq = qnorm(0.5 * (1 + conf[i]/100))
            lower[, i] = pred$pred - qq * pred$se
            upper[, i] = pred$pred + qq * pred$se}
        colnames(lower) = colnames(upper) = paste(conf, "%", sep = "")

        # Colors:
        shadecols = switch(1 + (length(conf) > 1), 7, length(conf):1)
        shadepalette = heat.colors(length(conf))
        col = 1

        # Plot History:
        npred = length(pred$pred)
        ylim = range(c(data[(n-n.back+1):n], pred$pred), na.rm = TRUE)
        ylim = range(ylim, lower, upper, na.rm = TRUE)
        ylab = paste("Series: ", fit$series)
        vTS = ts(c(data[(n-n.back+1):n], pred$pred[1], rep(NA, npred-1)),
            end = tsp(data)[2] + npred/freq, frequency = freq)
        plot(vTS, type = "o", pch = 19, ylim = ylim, ylab = ylab)
        title(main = paste(fit$tstitle))

        # Confidence Intervals:
        xx = tsp(data)[2] + (1:npred)/freq
        idx = rev(order(conf))
        if (nint > 1) palette(shadepalette)
        for (i in 1:nint) { polygon(c(xx, rev(xx)), c(lower[, idx[i]],
            rev(upper[, idx[i]])), col = shadecols[i], border = FALSE) }
        palette("default")

        # Mean:
        vTS = ts(pred$pred, start = tsp(data)[2]+1/freq, frequency = freq)
        lines(vTS, lty = 1, col = 4)
        points(vTS, pch = 19)

        # Printout:
        nconf = length(conf)
        out = pred$pred
        upper = as.matrix(upper)
        lower = as.matrix(lower)
        names = "Forecast"
        for (i in nconf:1) {
            out = cbind(out, lower[, i])
            names = c(names, paste("Low", conf[i])) }
        out = cbind(out, pred$pred)
        names = c(names, "Forecast")
        for (i in 1:nconf) {
            out = cbind(out, upper[, i])
            names = c(names, paste("High", conf[i])) }
        out = round(out, digits = 4)[,2:(2*nconf+2)]
        colnames(out) = names[2:(2*nconf+2)]

        # Grid:
        grid()
        options(warn = 0)

        # Add to Output:
        ans$out = out
    }

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.arPredict =
function (object, n.ahead = 10, se.fit = TRUE, ...)
{   # A function implemented by Diethelm Wuertz

    # FUNCTION:

    # Predict - object@fit$tsmodel = "ar":
    fit = object@fit
    class(fit) = "ar"
    ans = predict(object = fit, newdata = fit$x,
        n.ahead = n.ahead, se.fit = se.fit)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.arimaPredict =
function (object, n.ahead = 10, se.fit = TRUE, ...)
{   # A function implemented by Diethelm Wuertz

    # FUNCTION:

    # Predict - object@fit$tsmodel = "arima":
    fit = object@fit
    class(fit) = "Arima"
    if (!exists("xreg")) xreg = NULL
    if (!exists("newxreg")) newxreg = NULL
    class(object) = "Arima"
    ans = predict(object = fit, n.ahead = n.ahead,
        newxreg = newxreg, se.fit = se.fit, xreg = xreg, ...)

    # Return Value:
    ans
}


################################################################################

