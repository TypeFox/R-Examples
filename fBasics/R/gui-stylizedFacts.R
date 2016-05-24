
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
# You should have received A copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


###############################################################################
# FUNCTION:             STYLIZED FACTS GUI:
#  .stylizedFactsGUI     Opens a GUI for stylized facts
#  .lmacfPlot            Estimates and displays the long memory ACF
#  .logpdfPlot           Displays a pdf plot on logarithmic scale(s)
#  .ccfPlot              Displays tailored cross correlation function plot
#  .qqgaussPlot          Displays a tailored Gaussian quantile-quantile plot
###############################################################################


.lmacfPlot <-
function(x, lag.max = max(2, floor(10*log10(length(x)))),
    ci = 0.95, type = c("both", "acf", "hurst"), labels = TRUE,
    trace = TRUE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Evaluate and display long memory autocorrelation Function.

    # Arguments:
    #   x - an uni- or multivariate return series of class 'timeSeries'
    #       or any other object which can be transformed by the function
    #       'as.timeSeries()' into an object of class 'timeSeries'.
    #   labels - a logical flag, by default true. Should a default
    #       main title and labels addet to the plot?

    # FUNCTION:

    # Settings:
    if (!is.timeSeries(x)) x = as.timeSeries(x)
    Units = colnames(x)
    X = x

    # Match Arguments:
    type = match.arg(type)

    # Labels:
    if (labels) {
        main1 = ""
        xlab1 = "lag"
        ylab1 = "ACF"
        main2 = ""
        xlab2 = "log lag"
        ylab2 = "log ACF"
    } else {
        main1 = xlab1 = ylab1 = ""
        main2 = xlab2 = ylab2 = ""
    }

    Fitted = list()
    Hurst = NULL
    DIM = dim(X)[2]
    for (i in 1:DIM) {

        # Get Data:
        x.ret = as.matrix(X)[, i]
        x = abs(x.ret)
        if (labels) main1 = main2 = Units[i]

        # ACF:
        z = acf(x, lag.max = lag.max, type = "correlation", plot = FALSE)
        z$acf[1] = 0
        cl = qnorm(0.5 + ci/2)/sqrt(z$n.used)
        z.min = min(z$acf, -cl)

        # lin-lin plot excluding one:
        x = seq(0, lag.max, by = 1)
        y = z$acf
        if (type == "both" | type == "acf") {
            plot(x = x[-1], y = y[-1], type = "l", main = main1,
                col = "steelblue", xlab = xlab1, ylab = ylab1,
                xlim = c(0, lag.max), ylim = c(-2*cl, max(y[-1])), ...)
            # abline(h = 0, lty = 3)
            if (trace) {
                cat ("\nLong Memory Autocorrelation Function:")
                    paste (cat ("\n  Maximum Lag        "), cat(lag.max))
                    paste (cat ("\n  Cut-Off ConfLevel  "), cat(cl))
            }
            ACF = acf(x.ret, lag.max = lag.max, plot = FALSE)$acf[,,1]
            lines(x = 1:lag.max, y = ACF[-1], type = "l", col = "steelblue")
            lines(x = c(-0.1, 1.1)*lag.max, y = c(+cl, +cl), lty = 3,
                col = "darkgrey")
            lines(x = c(-0.1, 1.1)*lag.max, y = c(-cl, -cl), lty = 3,
                col = "darkgrey")
        }

        # log-log Plot of ACF:
        x = x[y > cl]
        y = y[y > cl]
        # log-log:
        if (length(x) < 10) {
            Fit = c(NA, NA)
            hurst = NA
            cat("\n  The time series exhibits no long memory! \n")
        } else {
            Fit = lsfit(log(x), log(y))
            fit = unlist(Fit)[1:2]
            hurst = 1 + fit[2]/2
            if (type == "both" | type == "hurst") {
                plot(x = log(x), y = log(y), type = "l", xlab = xlab2,
                    ylab = ylab2, main = main2, col = "steelblue", ...)
                Text = paste("Hurst Exponent:", signif(hurst, 3))
                mtext(Text, side = 4, adj = 0, col = "darkgrey", cex = 0.7)
                # Grid:
                if (labels) grid()
            }
            ### fit = l1fit(log(x), log(y))$coefficients
            abline(fit[1], fit[2], col = 1)
            if (trace) {
                paste (cat ('\n  Plot-Intercept     '), cat(fit[1]))
                paste (cat ('\n  Plot-Slope         '), cat(fit[2]))
                paste (cat ('\n  Hurst Exponent     '), cat(hurst), cat("\n"))
            }
        }

        # Save Results:
        if (DIM == 1) Fitted = Fit else Fitted[[i]] = Fit
        Hurst = c(Hurst, hurst)
    }

    # Return Value:
    invisible(list(fit = Fitted, hurst = Hurst))
}


# ------------------------------------------------------------------------------


.logpdfPlot =
function(x, breaks = "FD", type = c("lin-log", "log-log"),
doplot = TRUE, labels = TRUE, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays a pdf plot on logarithmic scale(s)

    # Details:
    #   Returns a pdf plot on a lin-log scale in
    #   comparison to a Gaussian density plot
    #   and return the break-midpoints and the
    #   counts obtained from the histogram of
    #   the empirical data.

    # Arguments:
    #   x - an uni- or multivariate return series of class 'timeSeries'
    #       or any other object which can be transformed by the function
    #       'as.timeSeries()' into an object of class 'timeSeries'.
    #   labels - a logical flag, by default true. Should a default
    #       main title and labels addet to the plot?

    # FUNCTION:

    # Settings:
    if (!is.timeSeries(x)) x = as.timeSeries(x)
    Units = colnames(x)

    # Select Type:
    type = match.arg(type)

    # Labels:
    if (labels) {
        if (type == "lin-log") {
            main = "log PDF"
            xlab = "x"
            ylab = "log PDF"
        } else if (type == "log-log") {
            main = "log PDF"
            xlab = "log x"
            ylab = "log PDF"
        }
    } else {
        main = xlab = ylab = ""
    }

    X = x
    DIM = ncol(X)

    for (i in 1:DIM) {

        # Get Data:
        x = as.vector(X[, i])
        if (labels) main = Units[i]

        # Lin-Log Plot:
        if (type == "lin-log") {

            # Histogram PDF:
            result = hist(x, breaks = breaks, plot = FALSE)
            prob.counts = result$counts/sum(result$counts) /
                diff(result$breaks)[1]
            histogram = list(breaks = result$breaks, counts = prob.counts)

            # Histogram Count & Break-Midpoints:
            yh = histogram$counts
            xh = histogram$breaks
            xh = xh[1:(length(xh)-1)] + diff(xh)/2
            xh = xh[yh > 0]
            yh = log(yh[yh > 0])
            if (doplot) {
                par(err = -1)
                plot(xh, yh, type = "p", pch = 19, col = "steelblue",
                    main = main, xlab = xlab, ylab = ylab, ...)
                Text = "Scales: log-log"
                mtext(Text, side = 4, adj =0, col = "darkgrey", cex = 0.7)
            }

            # Compare with a Gaussian Plot:
            xg = seq(from = xh[1], to = xh[length(xh)], length = 301)
            yg = log(dnorm(xg, mean(x), sqrt(var(x))))
            if (doplot) {
                par(err = -1)
                lines(xg, yg, col = "brown")
            }

            # Return Value:
            result = list(breaks = xh, counts = yh, fbreaks = xg,
                fcounts = yg)
        }

        # Log-Log Plot:
        if (type == "log-log") {

            # Histogram PDF:
            result = hist(x, breaks = breaks, plot = FALSE)
            prob.counts = result$counts/sum(result$counts) / diff(result$breaks)[1]
            histogram = list(breaks = result$breaks, counts = prob.counts)

            # Histogram Count & Breaks:
            yh = histogram$counts
            xh = histogram$breaks
            xh = xh[1:(length(xh)-1)] + diff(xh)/2
            xh = xh[yh > 0]
            yh = yh[yh > 0]
            yh1 = yh[xh < 0]
            xh1 = abs(xh[xh < 0])
            yh2 = yh[xh > 0]
            xh2 = xh[xh > 0]
            if (doplot) {
                plot(log(xh1), log(yh1), type = "p", pch = 19,
                    col = "darkgreen",
                    main = main, xlab = xlab, ylab = ylab, ...)
                Text = "Scales: log-log"
                mtext(Text, side = 4, adj =0, col = "darkgrey", cex = 0.7)
                par(err = -1)
                points(log(xh2), log(yh2), col = 2, ...)
            }

            # Compare with a Gaussian plot:
            xg = seq(from = min(xh1[1], xh[2]),
                to = max(xh1[length(xh1)], xh2[length(xh2)]), length = 301)
            xg = xg[xg > 0]
            yg = log(dnorm(xg, mean(x), sqrt(var(x))))
            if (doplot) {
                par(err = -1)
                lines(log(xg), yg, col = "brown")
            }

            # Result:
            result = list(breaks = c(xh1, xh2), counts = c(yh1, yh2),
                fbreaks = c(-rev(xg), xg), fcounts = c(-rev(yg), yg))
        }

        # Grid:
        if (labels) grid()

    }

    # Return Value:
    invisible(result)
}


# ------------------------------------------------------------------------------


.qqgaussPlot <-
function(x, span = 5, col = "steelblue", labels = TRUE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns a simple Quantile-Quantile plot.

    # Arguments:
    #   x - an uni- or multivariate return series of class 'timeSeries'
    #       or any other object which can be transformed by the function
    #       'as.timeSeries()' into an object of class 'timeSeries'.
    #   labels - a logical flag, by default true. Should a default
    #       main title and labels addet to the plot?

    # FUNCTION:

    # Settings:
    # if (SPLUSLIKE) splusLikePlot(TRUE)

    # Settings:
    if (!is.timeSeries(x)) x = as.timeSeries(x)
    Units = colnames(x)

    # Labels:
    if (labels) {
        main = "Normal QQ Plot"
        xlab = "Theoretical Quantiles"
        ylab = "Sample Quantiles"
    } else {
        main = xlab = ylab = ""
    }

    X = x
    DIM = dim(X)[2]

    for (i in 1:DIM) {

        # Get Data:
        x = as.vector(as.matrix(X)[, i])
        if (labels) main = Units[i]

        # Standardized qqnorm():
        y = (x-mean(x)) / sqrt(var(x))

        # Further Settings:
        y[abs(y) < span]
        lim = c(-span, span)

        # Plot qqnorm:
        qqnorm(y, main = main, xlab = xlab, ylab = ylab,
            xlim = lim, ylim = lim, col = col, ...)

        # Add Line:
        qqline(y, ...)

        # Grid:
        if (labels) grid()

    }

    # Return Value:
    invisible(x)
}


# ------------------------------------------------------------------------------


.ccfPlot <-
function(x, y, lag.max = max(2, floor(10*log10(length(x)))),
    type = c("correlation", "covariance", "partial"), labels = TRUE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Cross correlation function plot

    # Arguments:
    #   x - an univariate 'timeSeries' object
    #   labels - a logical flag, by default true. Should a default
    #       main title and labels addet to the plot?

    # FUNCTION:

    # Convert Type:
    if (class(x) == "timeSeries") stopifnot(isUnivariate(x))
    if (class(y) == "timeSeries") stopifnot(isUnivariate(y))
    x = as.vector(x)
    y = as.vector(y)

    # Labels:
    if (labels) {
        main = "Crosscorrelation Function"
        xlab = "lag"
        ylab = "CCF"
    } else {
        main = xlab = ylab = ""
    }

    # Result:
    # A copy from R's ccf - So you can use it also under SPlus:
    X = cbind(x = x, y = y)
    acf.out =  acf(X, lag.max = lag.max, plot = FALSE, type = type[1])
    lag = c(rev(acf.out$lag[-1, 2, 1]), acf.out$lag[, 1, 2])
    y = c(rev(acf.out$acf[-1, 2, 1]), acf.out$acf[, 1, 2])
    acf.out$acf = array(y, dim = c(length(y), 1, 1))
    acf.out$lag = array(lag, dim = c(length(y), 1, 1))
    acf.out$snames = paste(acf.out$snames, collapse = " & ")
    plot(acf.out, main = main, xlab = xlab, ylab = ylab, ...)

    # Return Value:
    invisible(acf.out)
}


################################################################################


.stylizedFactsGUI <-
function(x, mfrow = c(3, 3))
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Opens a GUI for stylized facts

    # FUNCTION:

    # Refresh Code:
    stylizedFactsRefreshCode <-
function(...)
    {
        # Settings:
        selectedAsset  = .tdSliderMenu(no = 1)
        type = as.integer(.tdSliderMenu(obj.name = "stylizedFactsType"))
        Unit = colnames(x)

        # ACF Plot:
        if (type == 1) {
            if (selectedAsset == 0) {
                par(mfrow = mfrow)
                acfPlot(x)
            } else {
                par(mfrow = c(1, 1))
                acfPlot(x[, selectedAsset])
            }
        }

        # PACF Plot:
        if (type == 2) {
            if (selectedAsset == 0) {
                par(mfrow = mfrow)
                pacfPlot(x)
            } else {
                par(mfrow = c(1, 1))
                pacfPlot(x[, selectedAsset])
            }
        }

        # Volatility ACF Plot:
        if (type == 3) {
            if (selectedAsset == 0) {
                par(mfrow = mfrow)
                acfPlot(abs(x))
            } else {
                par(mfrow = c(1, 1))
                acfPlot(abs(x[, selectedAsset]))
            }
        }

        # Taylor Effect Plot:
        if (type == 4) {
            if (selectedAsset == 0) {
                par(mfrow = mfrow)
                teffectPlot(x)
            } else {
                par(mfrow = c(1, 1))
                teffectPlot(x[, selectedAsset])
            }
        }

        # Long Memory ACF:
        if (type == 5) {
            if (selectedAsset == 0) {
                par(mfrow = mfrow)
                .lmacfPlot(abs(x))
            } else {
                par(mfrow = c(1, 1))
                .lmacfPlot(abs(x[, selectedAsset]))
            }
        }

        # Lagged Autocorrelation Plot:
        if (type == 6) {
            if (selectedAsset == 0) {
                par(mfrow = mfrow)
                lacfPlot(x)
            } else {
                par(mfrow = c(1, 1))
                lacfPlot(x[, selectedAsset])
            }
        }

        # PDF plot on lin-log Scale:
        if (type == 7) {
            if (selectedAsset == 0) {
                par(mfrow = mfrow)
                .logpdfPlot(x)
            } else {
                par(mfrow = c(1, 1))
                .logpdfPlot(x[, selectedAsset])
            }
        }

        # PDF plot on log-log Scale:
        if (type == 8) {
            if (selectedAsset == 0) {
                par(mfrow = mfrow)
                .logpdfPlot(x, type = "log-log")
            } else {
                par(mfrow = c(1, 1))
                .logpdfPlot(x[, selectedAsset], type = "log-log")
            }
        }

        # Simple Normal Quantile-Quantile Plot
        if (type == 9) {
            if (selectedAsset == 0) {
                par(mfrow = mfrow)
                .qqgaussPlot(x, pch = 19)
            } else {
                par(mfrow = c(1, 1))
                .qqgaussPlot(x[, selectedAsset], pch = 19)
            }
        }

        # Scaling Law Plot:
        if (type == 10) {
            if (selectedAsset == 0) {
                par(mfrow = mfrow)
                scalinglawPlot(x, pch = 19)
            } else {
                par(mfrow = c(1, 1))
                scalinglawPlot(x[, selectedAsset], pch = 19)
            }
        }

    }

    nAssets = dim(x)[2]

    .tdSliderMenu(
        stylizedFactsRefreshCode,

        names       = c("Selected Asset"),
        minima      = c(      0),
        maxima      = c(      nAssets),
        resolutions = c(      1),
        starts      = c(      0),

        but.functions = list(
        function(...){
                .tdSliderMenu(obj.name = "stylizedFactsType", obj.value = "1")
                stylizedFactsRefreshCode()},
        function(...){
                .tdSliderMenu(obj.name = "stylizedFactsType", obj.value = "2")
                stylizedFactsRefreshCode()},
        function(...){
                .tdSliderMenu(obj.name = "stylizedFactsType", obj.value = "3")
                stylizedFactsRefreshCode()},
        function(...){
                .tdSliderMenu(obj.name = "stylizedFactsType", obj.value = "4")
                stylizedFactsRefreshCode()},
        function(...){
                .tdSliderMenu(obj.name = "stylizedFactsType", obj.value = "5")
                stylizedFactsRefreshCode()},
        function(...){
                .tdSliderMenu(obj.name = "stylizedFactsType", obj.value = "6")
                stylizedFactsRefreshCode()},
        function(...){
                .tdSliderMenu(obj.name = "stylizedFactsType", obj.value = "7")
                stylizedFactsRefreshCode()},
        function(...){
                .tdSliderMenu(obj.name = "stylizedFactsType", obj.value = "8")
                stylizedFactsRefreshCode()},
        function(...){
                .tdSliderMenu(obj.name = "stylizedFactsType", obj.value = "9")
                stylizedFactsRefreshCode()},
        function(...){
                .tdSliderMenu(obj.name = "stylizedFactsType", obj.value = "10")
                stylizedFactsRefreshCode()}
        ),

        but.names = c(
            "1 ACF Function of Returns",
            "2 Partial ACF of Returns",
            "3 ACF of absolute Returns",
            "4 Taylor Effect",
            "5 Long Memory ACF of abs Returns",
            "6 Lagged Autocorrelations",
            "7 Lin-Log Tail Density Plot",
            "8 Log-Log Lower Tail Density",
            "9 Simple Normal QQ Plot",
            "10 Scaling Law Plot"),

        title = "Stylized Facts GUI"
        )

   .tdSliderMenu(obj.name = "type", obj.value = "1", no = 1)

   # return Value:
   invisible()
}


################################################################################

