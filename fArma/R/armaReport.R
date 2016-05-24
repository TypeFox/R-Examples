
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
#  show.fARMA              S4: Prints a fitted ARMA time series object
#  plot.fARMA              S3: Plots stylized facts of a fitted ARMA object
#  summary.fARMA           S3: Summarizes a fitted ARMA time series object
################################################################################


################################################################################
# BUILTIN - PACKAGE DESCRIPTION:
#  Package: fracdiff
#  Version: 1.1-1
#  Title: Fractionally differenced ARIMA (p,d,q) models
#  Date: 2004-01-12
#  Author: S original by Chris Fraley <fraley@stat.washington.edu>.
#    R port by Fritz Leisch <leisch@ci.tu-wien.ac.at>;
#    since 2003-12: Martin Maechler
#  Maintainer: Martin Maechler <maechler@stat.math.ethz.ch>
#  Description: Maximum likelihood estimation of the parameters of a
#    fractionally differenced ARIMA(p,d,q) model (Haslett and Raftery,
#    Appl.Statistics, 1989).
#  License: GPL version 2 or later
#  Packaged: Mon Jan 12 11:22:27 2004; maechler
################################################################################


setClass("fARMA",
    representation(
        call = "call",
        formula = "formula",
        method = "character",
        parameter = "list",
        data = "list",
        fit = "list",
        residuals = "list",
        fitted = "list",
        title = "character",
        description = "character"
    )
)


# ------------------------------------------------------------------------------


setMethod("show", "fARMA",
    function(object)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Old S3 print method for a fitted ARMA timeSeries object

    # FUNCTION:

    # Title:
    cat("\nTitle:\n ")
    cat(object@title, "\n")

    # Call:
    cat("\nCall:\n ")
    cat(paste(deparse(object@call), sep = "\n", collapse = "\n"),
        "\n", sep = "")

    # Model:
    cat("\nModel:\n ", object@fit$tstitle, "\n", sep = "")

    # Coefficients:
    cat("\nCoefficient(s):\n")
    digits = max(4, getOption("digits") - 4)
    print.default(format(object@fit$coef, digits = digits), print.gap = 2,
        quote = FALSE)

    # Description:
    cat("\nDescription:\n ")
    cat(object@description, "\n\n")

    # Return Value:
    invisible()
})


# ------------------------------------------------------------------------------





# ------------------------------------------------------------------------------


summary.fARMA =
function (object, doplot = TRUE, which = "all", ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Analyzes a Fitted ARMA timeSeries Object

    # FUNCTION:

    # Initialize:
    if (object@fit$tsmodel == "arfima" & doplot) {
        warning(" Plot Method for arfima Models not yet implemented")
        doplot = FALSE
    }
    ans = NULL

    # Fit Call and Model:
    x = object
    object = x@fit
    ans$call = object$call
    ans$tsmodel = object$tstitle

    # Calculate Residuals and Variance:
    # ans$residuals = na.remove(object$residuals)
    ans$residuals = as.vector(na.omit(object$residuals))
    if (length(ans$residuals) == 0) {
        ans$var = 0 }
    if (length(ans$residuals) > 0) {
        ans$var = var(ans$residuals) }
    ans$sigma2 = object$sigma2

    # Generate Coef Matrix:
    tval = object$coef/object$se.coef
    prob = 2 * (1 - pnorm(abs(tval)))
    ans$coefmat = cbind(object$coef, object$se.coef, tval, prob)
    dimnames(ans$coefmat) = list(names(object$coef),
        c(" Estimate", " Std. Error", " t value", "Pr(>|t|)"))

    # More Parameters: aic, etc ...
    if (object$tsmodel == "ar") {
        ans$aic = (object$n.used * (1 + log(2 * pi)) + object$n.used *
            log(ans$var) + 2 * length(object$coef)) }
    if (object$tsmodel == "arma") {
        ans$aic = (object$n.used * (1 + log(2 * pi)) + object$n.used *
            log(ans$var) + 2 * length(object$coef))
        ans$css = object$css }
    if (object$tsmodel == "arima") {
        ans$aic = object$aic
        ans$loglik = object$loglik }
    if (object$tsmodel == "fracdiff") {
        doplot = FALSE }

    # Print Part:

    # Title:
    cat("\nTitle:\n ")
    cat(x@title, "\n")

    # Call:
    cat("\nCall:\n ")
    cat(paste(deparse(object$call), sep = "\n", collapse = "\n"),
        "\n", sep = "")

    # Model:
    cat("\nModel:\n ", object$tstitle, "\n", sep = "")

    # Coefficients:
    cat("\nCoefficient(s):\n")
    digits = max(4, getOption("digits") - 4)
    print.default(format(object$coef, digits = digits), print.gap = 2,
        quote = FALSE)

    # Residuals:
    digits = max(4, getOption("digits") - 4)
    if (length(object$residuals) > 2) {
        cat("\nResiduals:\n")
        rq = structure(quantile(ans$residuals),
            names = c("Min", "1Q", "Median", "3Q", "Max"))
        print(rq, digits = digits)
        # Moments:
        cat("\nMoments: \n")
        skewness = sum((ans$residuals - mean(ans$residuals))^3 /
            sqrt(var(ans$residuals))^3)/length(ans$residuals)
        kurtosis = sum((ans$residuals - mean(ans$residuals))^4 /
            var(ans$residuals)^2)/length(ans$residuals) - 3
        stats = structure(c(skewness, kurtosis),
            names = c("Skewness", "Kurtosis"))
        print(stats, digits = digits) }

    # Coef Matrix:
    cat("\nCoefficient(s):\n")
    signif.stars = getOption("show.signif.stars")
    printCoefmat(ans$coefmat, digits = digits,
        signif.stars = signif.stars, ...)

    # Fit:
    cat("\n")
    if (x@fit$tsmodel == "ar") {
        cat("sigma^2 estimated as:       ",
            format(object$var, digits = digits), "\n")
        cat("AIC Criterion:              ",
            format(round(object$aic, 2)), "\n") }
    if (x@fit$tsmodel == "arma") {
        cat("sigma^2 estimated as:       ",
            format(object$sigma2, digits = digits), "\n")
        cat("Conditional Sum-of-Squares: ",
            format(round(object$css, digits=2)), "\n")
        ## cat("AIC Criterion:              ",
        ##    format(round(object$aic, digits=2)), "\n")
        }
    if (x@fit$tsmodel == "arima") {
        cm = object$call$method
        if (is.null(cm) || cm != "CSS")
            cat(
              "sigma^2 estimated as: ", format(object$sigma2, digits = digits),
            "\nlog likelihood:       ", format(round(object$loglik, 2)),
            "\nAIC Criterion:        ", format(round(object$aic, 2)),
            "\n", sep = "")
        else
            cat(
              "sigma^2 estimated as: ", format(object$sigma2, digits = digits),
            "\npart log likelihood:  ", format(round(object$loglik,2)),
            "\n", sep = "") }

    # Doplot:
    if (doplot) plot.fARMA(x, which = which, ...)

    # Description:
    cat("\nDescription:\n ")
    cat(x@description, "\n\n")

    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


plot.fARMA =
function(x, which = "ask", gof.lag = 10, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Plot method for an object of class 'fARMA'

    # FUNCTION:

    # Check:
    if (x@fit$tsmodel == "arfima") {
        warning(" Plot method for ARFIMA models not yet implemented")
        return()
    }

    # Store Lag:
    x@fit$gof.lag = gof.lag

    # Plot:
    .interactiveArmaPlot(
        x,
        choices = c(
            "Standardized Residuals",
            "ACF of Residuals",
            "QQ Plot of Residuals",
            "Ljung-Box p Values"),
        plotFUN = c(
            ".plot.arma.1",  ".plot.arma.2",  ".plot.arma.3", ".plot.arma.4"),
        which = which)

    # Return Value:
    invisible(x)
}


# ------------------------------------------------------------------------------


.plot.arma.1 <-
function(x, ...)
{
    # 1. Standardized Residuals Plot:
    object = x@fit
    rs = as.vector(na.omit(object$residuals))
    stdres = rs/sqrt(object$sigma2)
    plot(stdres, type = "h",
        main = "Standardized Residuals",
        ylab = "Residuals", col = "steelblue", ...)
    grid()
    abline(h = 0, col = "grey")
}


# ------------------------------------------------------------------------------


.plot.arma.2 <-
function(x, ...)
{
    # 2. ACF of Residuals:
    object = x@fit
    acf(object$residuals, plot = TRUE, main = "ACF of Residuals",
        na.action = na.pass, ...)
    grid()
}


# ------------------------------------------------------------------------------


.plot.arma.3 <-
function(x, ...)
{
    # 3. QQ Plot of Residuals:
    object = x@fit
    rs = as.vector(na.omit(object$residuals))
    stdres = rs/sqrt(object$sigma2)
    qqnorm(stdres,
        xlab = "Normal Quantiles",
        ylab = "Residual Quantiles",
        main = "QQ-Plot of Residuals",
        pch = 19, col = "steelblue", ...)
    qqline(stdres, col = "grey")
    grid()
}


# ------------------------------------------------------------------------------


.plot.arma.4 <-
function(x, ...)
{
    # 4. Ljung-Box p Values:
    object = x@fit
    rs = as.vector(na.omit(object$residuals))
    nlag = x@fit$gof.lag
    pval = numeric(nlag)
    for (i in 1:nlag)
        pval[i] = Box.test(rs, i, type = "Ljung-Box")$p.value
    plot(1:nlag, pval, xlab = "lag", ylab = "p value", ylim = c(0, 1),
        pch = 19, col = "steelblue", main = "Ljung-Box p-values", ...)
    abline(h = 0.05, lty = 2, col = "grey")
    grid()
}


# ------------------------------------------------------------------------------

.interactiveArmaPlot =
function(x, choices = paste("Plot", 1:19),
plotFUN = paste("plot.", 1:19, sep = ""), which = "all", ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Plot method for an object of class "template".

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

    # Note:
    #   At maximum 19 plots are supported.

    # FUNCTION:

    # Some cecks:
    if (length(choices) != length(plotFUN))
        stop("Arguments choices and plotFUN must be of same length.")
    if (length(which) > length(choices))
        stop("Arguments which has incorrect length.")
    if (length(which) > length(plotFUN))
        stop("Arguments which has incorrect length.")
    if (length(choices) > 19)
        stop("Sorry, only 19 plots at max are supported.")

    # Plot:
    if (is.numeric(which)) {
        Which = rep(FALSE, times = length(choices))
        Which[which] = TRUE
        which = Which
    }
    if (which[1] == "all") {
        which = rep(TRUE, times = length(choices))
    }
    if (which[1] == "ask") {
        .multArmaPlot(x, choices, ...)
    } else {
        for ( i in 1:length(which) ) {
            FUN = match.fun(plotFUN[i])
            if (which[i]) FUN(x)
        }
    }

    # Return Value:
    invisible(x)
}


# ------------------------------------------------------------------------------


.multArmaPlot =
function (x, choices, ...)
{
    # Internal "askPlot" Function:

    # Match Functions, up to nine ...
    if (length(plotFUN) < 19) plotFUN =
        c(plotFUN, rep(plotFUN[1], times = 19 - length(plotFUN)))
    plot.1  = match.fun(plotFUN[1]);  plot.2  = match.fun(plotFUN[2])
    plot.3  = match.fun(plotFUN[3]);  plot.4  = match.fun(plotFUN[4])
    plot.5  = match.fun(plotFUN[5]);  plot.6  = match.fun(plotFUN[6])
    plot.7  = match.fun(plotFUN[7]);  plot.8  = match.fun(plotFUN[8])
    plot.9  = match.fun(plotFUN[9]);  plot.10 = match.fun(plotFUN[10])
    plot.11 = match.fun(plotFUN[11]); plot.12 = match.fun(plotFUN[12])
    plot.13 = match.fun(plotFUN[13]); plot.14 = match.fun(plotFUN[14])
    plot.15 = match.fun(plotFUN[15]); plot.16 = match.fun(plotFUN[16])
    plot.17 = match.fun(plotFUN[17]); plot.18 = match.fun(plotFUN[18])
    plot.19 = match.fun(plotFUN[19])
    pick = 1
    while (pick > 0) { pick = menu (
        ### choices = paste("plot:", choices),
        choices = paste(" ", choices),
        title = "\nMake a plot selection (or 0 to exit):")
        # up to 19 plot functions ...
        switch (pick,
            plot.1(x),  plot.2(x),  plot.3(x),  plot.4(x),  plot.5(x),
            plot.6(x),  plot.7(x),  plot.8(x),  plot.9(x),  plot.10(x),
            plot.11(x), plot.12(x), plot.13(x), plot.14(x), plot.15(x),
            plot.16(x), plot.17(x), plot.18(x), plot.19(x))
    }
}


################################################################################

