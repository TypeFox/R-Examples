
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA

# Copyrights (C)
# for this R-port:
#   1999 - 2008, Diethelm Wuertz, Rmetrics Foundation, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:                 DESCRIPTION:
#  plot                      Plot method for an object of class 'fGARCH'
#  .interactiveGarchPlot     Plot interactively
#  .multGarchPlot            Arrange multivariate Plots
#  .plot.garch.1               Plot Time Series
#  .plot.garch.2               Plot Conditional SD
#  .plot.garch.3               Plot Series with 2 Conditional SD Superimposed
#  .plot.garch.4               Plot ACF of Observations
#  .plot.garch.5               Plot ACF of Squared Observations
#  .plot.garch.6               Plot Cross Correlation
#  .plot.garch.7               Plot Residuals
#  .plot.garch.8               Plot Conditional SDs
#  .plot.garch.9               Plot Standardized Residuals
#  .plot.garch.10              Plot ACF of Standardized Residuals
#  .plot.garch.11              Plot ACF of Squared Standardized Residuals
#  .plot.garch.12              Plot Cross Correlation between r^2 and r
#  .plot.garch.13              Plot QQ-Plot of Standardized Residuals"
#   .qqDist                     Quantile-Quantile Points
#   .qqLine                     Quantile-Quantile Line
################################################################################


setMethod(f = "plot", signature(x = "fGARCH", y = "missing"), definition =
    function(x, which = "ask", ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Plot method for an object of class 'fGARCH'

    # Note:
    #   This method can also be used for plotting graphs fitted by
    #   the function 'garch' from the contributed R package 'tseries'.

    # FUNCTION:

    if (as.character(x@call[1]) == ".gogarchFit") 
    {
        # Plot multivariate GO-Garch model:
        print("GO-Garch Plot Not Yet Implemented")
    } else {
        # Plot univariate Models:
        .interactiveGarchPlot(
            x,
            choices = c(
                "Time Series",
                "Conditional SD",
                "Series with 2 Conditional SD Superimposed",
                "ACF of Observations",
                "ACF of Squared Observations",
                "Cross Correlation",
                "Residuals",
                "Conditional SDs",
                "Standardized Residuals",
                "ACF of Standardized Residuals",
                "ACF of Squared Standardized Residuals",
                "Cross Correlation between r^2 and r",
                "QQ-Plot of Standardized Residuals"),
            plotFUN = paste(".plot.garch", 1:13, sep = "."),
            which = which, ...)
    }
    
    # Return Value:
    invisible(x)
})


# ------------------------------------------------------------------------------


.interactiveGarchPlot <-
    function(x, choices, plotFUN, which, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Internal plot function

    # Arguments:
    #   x - an object to be plotted
    #   choices - the character string for the choice menu
    #   plotFUN - the names of the plot functions
    #   which - plot selection, which graph should be
    #     displayed. If a character string named "ask" the
    #     user is interactively asked which to plot, if
    #     a logical vector of length N, those plots which
    #     are set "TRUE" are displayed, if a character string
    #     named "all" all plots are displayed.

    # FUNCTION:

    # Some cecks:
    if (length(choices) != length(plotFUN))
        stop("Arguments choices and plotFUN must be of same length.")
    if (length(which) > length(choices))
        stop("Arguments which has incorrect length.")
    if (length(which) > length(plotFUN))
        stop("Arguments which has incorrect length.")

    # Plot:
    if (is.numeric(which)) {
        Which = rep(FALSE, times = length(choices))
        Which[which] = TRUE
    }

    if (which[1] == "all") {
        Which = rep(TRUE, times = length(choices))
    }

    if (which[1] == "ask") {
        .multGarchPlot(x, choices, plotFUN, ...)
    } else {
        for ( i in 1:length(choices) ) {
            FUN = match.fun(plotFUN[i])
            if (Which[i]) FUN(x)
        }
    }

    # Return Value:
    invisible(x)
}


# ------------------------------------------------------------------------------


.multGarchPlot <-
    function (x, choices, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Internal plot function

    pick = 1
    while (pick > 0) {
        pick = menu (
            ### choices = paste("plot:", choices),
            choices = paste(" ", choices),
            title = "\nMake a plot selection (or 0 to exit):")
        # up to 19 plot functions ...
        switch (pick,
            .plot.garch.1(x),  .plot.garch.2(x),  .plot.garch.3(x),
            .plot.garch.4(x),  .plot.garch.5(x),  .plot.garch.6(x),
            .plot.garch.7(x),  .plot.garch.8(x),  .plot.garch.9(x),
            .plot.garch.10(x), .plot.garch.11(x), .plot.garch.12(x),
            .plot.garch.13(x))
    }
}


# ------------------------------------------------------------------------------


.plot.garch.1 <-
    function(x, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Internal plot function

    # 1. Time Series:
    xseries = x@data
    plot(xseries, type = "l", col = "steelblue", ylab = "x",
        main = "Time Series")
    abline(h = 0, col = "grey", lty = 3)
    grid()
}


# ------------------------------------------------------------------------------


.plot.garch.2 <-
    function(x, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Internal plot function

    # 2. Conditional SD:
    xcsd = volatility(x, "sigma")
    plot(xcsd, type = "l", col = "steelblue", ylab = "x",
        main = "Conditional SD")
    abline(h = 0, col = "grey", lty = 3)
    grid()
}


# ------------------------------------------------------------------------------


.plot.garch.3 <-
    function(x, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Internal plot function

    # 3. Series with 2 Conditional SD Superimposed:
    xseries = x@data
    xcsd = volatility(x, "sigma")
    ci = 2
    plot(xseries, type = "l", col = "steelblue", ylab = "x",
        main = "Series with 2 Conditional SD Superimposed")
    lines(mean(xseries) + ci * xcsd, col = "grey") # or simply xseries ?
    lines(mean(xseries) - ci * xcsd, col = "grey")
    abline(h = 0, col = "grey", lty = 3)
    grid()
}


# ------------------------------------------------------------------------------


.plot.garch.4 <-
    function(x, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Internal plot function

    # 4. ACF of the Observations:
    xseries = as.vector(x@data)
    n = length(xseries)
    lag.max = as.integer(10*log10(n))
    acf(xseries, lag.max = lag.max, xlab = "Lags", col = "steelblue",
        main = "ACF of Observations", plot = TRUE)
}


# ------------------------------------------------------------------------------


.plot.garch.5 <-
    function(x, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Internal plot function

    # 5. ACF of the Squared Observations:
    xseries = as.vector(x@data)
    xseries2 = xseries^2
    n = length(xseries)
    lag.max = as.integer(10*log10(n))
    acf(xseries2, lag.max = lag.max, xlab = "Lags", col = "steelblue",
        main = "ACF of Squared Observations", plot = TRUE)
}


# ------------------------------------------------------------------------------


.plot.garch.6 <-
    function(x, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Internal plot function

    # 6. Cross Correlation between x^2 and x:
    xseries = as.vector(x@data)
    xseries2 = xseries^2
    n = length(xseries)
    lag.max = as.integer(10*log10(n))
    ccf(xseries2, xseries, lag.max = lag.max, xlab = "Lags",
        main = "Cross Correlation", plot = TRUE, col = "steelblue")
}


# ------------------------------------------------------------------------------


.plot.garch.7 <-
    function(x, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Internal plot function

    # 7. Residuals:
    res = residuals(x, standardize = FALSE)
    plot(res, type = "l", main = "Residuals", col = "steelblue", ...)
    abline(h = 0, lty = 3)
    grid()
}


# ------------------------------------------------------------------------------


.plot.garch.8 <-
    function(x, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Internal plot function

    # 8. Conditional SDs:
    xcsd = volatility(x, "sigma")
    plot(xcsd, type = "l", main = "Conditional SD's",
        col = "steelblue", ...)
    abline(h = 0, lty = 3)
    grid()
}


# ------------------------------------------------------------------------------


.plot.garch.9 <-
    function(x, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Internal plot function

    # 9. Standardized Residuals:
    sres = residuals(x, standardize = TRUE)
    plot(sres, type = "l", main = "Standardized Residuals",
        col = "steelblue", ...)
    abline(h = 0, lty = 3)
    grid()
}


# ------------------------------------------------------------------------------


.plot.garch.10 <-
    function(x, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Internal plot function

    # 10. ACF of Standardized Residuals:
    sres = as.matrix(residuals(x, standardize = TRUE))
    n = length(sres)
    lag.max = as.integer(10*log10(n))
    acf(sres, lag.max = lag.max, xlab = "Lags", col = "steelblue",
        main = "ACF of Standardized Residuals", plot = TRUE)
}


# ------------------------------------------------------------------------------


.plot.garch.11 <-
    function(x, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Internal plot function

    # 11. ACF of Squared Standardized Residuals:
    sres2 = as.matrix(residuals(x, standardize = TRUE)^2)
    n = length(sres2)
    lag.max = as.integer(10*log10(n))
    acf(sres2, lag.max = lag.max, xlab = "Lags", col = "steelblue",
        main = "ACF of Squared Standardized Residuals", plot = TRUE)
}


# ------------------------------------------------------------------------------


.plot.garch.12 <-
function(x, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Internal plot function

    # 12. Cross Correlation between r^2 and r:
    sres = residuals(x, standardize = FALSE)
    sres2 = sres^2
    n = length(sres)
    lag.max = as.integer(10*log10(n))
    ccf(sres2, sres, lag.max = lag.max, xlab = "Lags",
        main = "Cross Correlation", plot = TRUE, col = "steelblue")
}


# ------------------------------------------------------------------------------


.plot.garch.13 <-
    function(x, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Internal plot function

    # 13. QQ-Plot of Standardized Residuals:
    sres = residuals(x, standardize = TRUE)
    cond.dist = x@fit$params$cond.dist
    cond.dist = paste("q", cond.dist, sep = "")
    nc = nchar(x@fit$params$cond.dist)

    parNames <- names(x@fit$par)
    skew <-
        if ("skew" %in% parNames)
            x@fit$par["skew"]
        else
            x@fit$params$skew
    shape <-
        if ("shape" %in% parNames)
            x@fit$par["shape"]
        else
            x@fit$params$shape

    if (cond.dist == "qnorm" || cond.dist == "qQMLE")
        .qqDist(sres, dist = "qnorm")

    if (cond.dist == "qstd" | cond.dist == "qged")
        .qqDist(sres, dist = cond.dist, nu = shape)

    if (cond.dist == "qsnorm")
        .qqDist(sres, dist = cond.dist, xi = skew)

    if (cond.dist == "qsstd" | cond.dist == "qsged")
        .qqDist(sres, dist = cond.dist, xi = skew, nu = shape)

    if (cond.dist == "qsnig")
        .qqDist(sres, dist = ".qsnigC", rho = skew, zeta = shape)

}


# ------------------------------------------------------------------------------


.qqDist <-
    function (y, dist = "qnorm", ylim = NULL, main = paste(dist, "- QQ Plot"),
    xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", doplot = TRUE,
    datax = FALSE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description
    #   QQ Plot for arbitray distribution

    # FUNCTION:
    # print(dist)

    # Match Function :
    qDist = match.fun(dist)

    # Check Arguments:
    # if (substr(dist, 1, 1) != "q") stop("dist is misspecified")
    # test = class(test = try(qDist(0.5, ...), silent = TRUE))
    # if (test == "try-error") stop("dist does not exist")

    # Transform to Vector Mode:
    y = as.vector(y)

    # Compute Data:
    if (has.na <- any(ina <- is.na(y))) {
        yN = y
        y = y[!ina]
    }
    if (0 == (n <- length(y))) stop("y is empty or has only NAs")
    x <- qDist(ppoints(n,), ...)[order(order(y))]
    if (has.na) {
        y = x
        x = yN
        x[!ina] = y
        y = yN
    }

    # Create QQ Plot:
    if (doplot) {
        if (is.null(ylim)) ylim = range(y)
        if (datax) {
            plot(y, x, main = main, xlab = ylab, ylab = xlab, xlim = ylim,
                col = "steelblue", cex = 0.7)
        } else {
            plot(x, y, main = main, xlab = xlab, ylab = ylab, ylim = ylim,
                col = "steelblue", cex = 0.7)
        }
        .qqLine(y = y, dist = dist, datax = datax, ...)
        grid()
    }

    # Return Value:
    invisible(if (datax) list(x = y, y = x) else list(x = x, y = y))
}


# ------------------------------------------------------------------------------


.qqLine <-
function (y, dist = "qnorm", datax = FALSE, ...)
{   # A function implemented by Diethelm Wuertz

    # Description
    #   Add slope to QQ Plot for arbitray distribution

    # FUNCTION:

    # Match Function :
    qDist = match.fun(dist)

    # Check Arguments:
    # if (substr(dist, 1, 1) != "q") stop("dist is misspecified")
    # test = class(test = try(qDist(0.5, ...), silent = TRUE))
    # if (test == "try-error") stop("dist does not exist")

    # Transform to Vector Mode:
    y = as.vector(y)

    # Compute Data:
    y = quantile(y[!is.na(y)], c(0.25, 0.75))
    x = qDist(c(0.25, 0.75), ...)

    # Add Slope:
    if (datax) {
        slope <- diff(x)/diff(y)
        int <- x[1] - slope * y[1]
    } else {
        slope <- diff(y)/diff(x)
        int <- y[1] - slope * x[1]
    }

    # Return Value:
    abline(int, slope)
}


################################################################################

