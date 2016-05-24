
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
# METHODS:              PRINT, PLOT, AND SUMMARY:
#  show.fGEVFIT          S4 Show method for object of class "fGEVFIT"
#  plot.fGEVFIT          S3 Plot method for object of class "fGEVFIT"
#   .gev1Plot             Block Maxima Plot
#   .gev2Plot             Scatterplot of Residuals
#   .gev3Plot             Histogram of Residuals
#   .gev4Plot             Quantile-Quantile Plot
#  summary.fGEVFIT       S3 Summary Method for object of class "fGEVFIT"
################################################################################


setMethod("show", "fGEVFIT",
    function(object)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Print Method for an object of class "fGEVFIT".

    # Arguments:
    #   object - an object of class fGEVFIT

    # FUNCTION:

    # Title:
    cat("\nTitle:\n" , object@title, "\n")

    # Function Call:
    cat("\nCall:\n ")
    cat(paste(deparse(object@call), sep = "\n",
        collapse = "\n"), "\n", sep = "")

    # Estimation Type:
    cat("\nEstimation Type:\n ", object@method, "\n")

    # Estimated Parameters:
    cat("\nEstimated Parameters:\n")
    print(object@fit$par.ests)

    # Desription:
    cat("\nDescription\n ", object@description, "\n\n")

    # Return Value:
    invisible(object)
})





################################################################################


plot.fGEVFIT =
    function(x, which = "ask", ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Plot method for an object of class "gevFit".

    # Arguments:

    # Details:
    #   plot.gev:
    #   Data are converted to unit exponentially distributed residuals
    #   under null hypothesis that GEV fits. Two diagnostics for iid
    #   exponential data are offered:
    #   "Scatterplot of Residuals" and "QQplot of Residuals"

    # Example:
    #   fit = gevFit(gevSim(), type = "mle", gumbel = FALSE)
    #   par(mfrow = c(2, 2)); plot(fit)
    #   par(mfrow = c(1, 1)); plot(fit, which = "ask")
    #
    #   fit = gevFit(gevSim(), type = "mle", gumbel = TRUE)
    #   par(mfrow = c(1, 1)); plot(fit, which = "ask")
    #
    #   fit = gevFit(gevSim(), type = "pwm", gumbel = FALSE)
    #   par(mfrow = c(1, 1)); plot(fit, which = "ask")
    #
    #   fit = gevFit(gevSim(), type = "pwm", gumbel = TRUE)
    #   par(mfrow = c(1, 1)); plot(fit, which = "ask")

    # FUNCTION:

    # Plot:
    interactivePlot(
        x = x,
        choices = c(
            "Block Maxima Plot",
            "Scatterplot of Residuals",
            "Histogram of Residuals",
            "Quantile Quantile Plot"),
        plotFUN = c(
            ".gev1Plot",
            ".gev2Plot",
            ".gev3Plot",
            ".gev4Plot"),
        which = which)

    # Return Value:
    invisible(x)
}


# ------------------------------------------------------------------------------


.gev1Plot =
    function(x, labels = TRUE, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Time Series Plot of Block Maxima

    # Arguments:
    #   x - an object of class "fGEVFIT" as returned by the 
    #       function gevFit
    #   labels - a logical, should labels be added to the plot
    #   ... - optional arguments passed to the function plot
    
    # Example:
    #   .gev1Plot(gevFit(gevSim()))

    # FUNCTION:

    # Data:
    data = x@data$blockmaxima

    # Labels:
    if (labels) {
        main = "Block Maxima"
        xlab = "Index"
        ylab = "Data"
    } else {
        main = xlab = ylab = ""
    }

    # Plot:
    plot(data, type = "h",
        main = main, xlab = xlab, ylab = ylab,
        col = "steelblue", ...)

    # Add Grid:
    if (labels) grid()

    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


.gev2Plot =
    function(x, labels = TRUE, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Scatterplot of Residuals:

    # Arguments:
    #   x - an object of class "fGEVFIT" as returned by the 
    #       function gevFit
    #   labels - a logical, should labels be added to the plot
    #   ... - optional arguments passed to the function plot
    
    # Example:
    #   .gev2Plot(gevFit(gevSim()))

    # FUNCTION:

    # Data:
    residuals = x@residuals

    # Labels:
    if (labels) {
        main = "Scatterplot of Residuals"
        xlab = "Ordering"
        ylab = "Residuals"
    } else {
        main = xlab = ylab = ""
    }

    # Plot:
    plot(residuals,
        main = main, xlab = xlab, ylab = ylab,
        pch = 19, col = "steelblue", ...)
    lines(lowess(1:length(residuals), residuals), col = "brown")

    # Add Grid:
    if (labels) grid()

    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


.gev3Plot =
    function(x, labels = TRUE, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Histogram Plot of Residuals with Gaussian Fit:

    # Arguments:
    #   x - an object of class "fGEVFIT" as returned by the 
    #       function gevFit
    #   labels - a logical, should labels be added to the plot
    #   ... - optional arguments passed to the function hist
    
    # Example:
    #   .gev3Plot(gevFit(gevSim()))

    # FUNCTION:

    # Data:
    residuals = x@residuals

    # Labels:
    if (labels) {
        if (x@method[1] == "gev") {
            dist = "GEV"
        } else if (x@method[1] == "gum") {
            dist = "Gumbel"
        }
        main = paste(dist, "Residual Histogram")
        xlab = "Residuals"
        ylab = "Density"
    } else {
        main = xlab = ylab = ""
    }

    # Plot:
    hist(residuals, probability = TRUE, breaks = "FD",
        main = main, xlab = xlab, ylab = ylab,
        col = "steelblue", border = "white", ...)

    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


.gev4Plot =
function(x, labels = TRUE, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Quantile-Quantile Plot:

    # Arguments:
    #   x - an object of class "fGEVFIT" as returned by the 
    #       function gevFit
    #   labels - a logical, should labels be added to the plot
    #   ... - optional arguments passed to the function plot
    
    # Example:
    #   .gev4Plot(gevFit(gevSim()))

    # FUNCTION:

    # Data:
    data = x@residuals
    sorted = sort(data)
    y <- qexp(ppoints(data))

    # Labels:
    if (labels) {
        main = "QQ Plot of Residuals"
        xlab = "Ordered Data"
        ylab = "Exponential Quantiles"
    } else {
        main = xlab = ylab = ""
    }

    # Plot:ata, type = "h",
    plot(sorted, y,
        main = main, xlab = xlab, ylab = ylab,
        pch = 19, col = "steelblue", ...)
    abline(lsfit(sorted, y))

    # Add Grid:
    if (labels) grid()

    # Return Value:
    invisible()
}


################################################################################


summary.fGEVFIT =
    function(object, doplot = TRUE, which = "all", ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Summary method for an object of class "gevFit".

    # Arguments:
    #   object - an object of class "fGEVFIT" as returned by the 
    #       function gevFit
    #   doplot - a logical, should a plot be returned 
    #   which - which plot(s) should be returned 
    #   optional arguments passed to the function plot
    
    # Example:
    #   fit = gevFit(gevSim(), type = "mle")
    #   par(mfrow = c(2, 2)); summary(fit)

    # FUNCTION:

    # Title:
    cat("\nTitle:\n", object@title, "\n")

    # Function Call:
    cat("\nCall:\n ")
    cat(paste(deparse(object@call), sep = "\n",
        collapse = "\n"), "\n", sep = "")

    # Estimation Type:
    cat("\nEstimation Type:\n ", object@method, "\n")

    # Estimated Parameters:
    cat("\nEstimated Parameters:\n")
    print(object@fit$par.ests)

    # Summary:
    if (object@method[2] == "mle") {
        cat("\nStandard Deviations:\n ");
        print(object@fit$par.ses)
        cat("\nLog-Likelihood Value:\n ", object@fit$llh, "\n")
        cat("\nType of Convergence:\n ", object@fit$converged, "\n") }

    # Desription:
    cat("\nDescription\n ", object@description, "\n\n")

    # Plot:
    if (doplot) {
        plot(object, which = which, ...)
    }

    # Return Value:
    invisible(object)
}


################################################################################

