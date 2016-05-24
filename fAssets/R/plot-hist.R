
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


################################################################################
# FUNCTION:                   DESCRIPTION:
#  assetsHistPlot              Displays a histograms of a single asset 
#  assetsLogDensityPlot        Displays a pdf plot on logarithmic scale
################################################################################


assetsHistPlot =
    function(x, col = "steelblue", skipZeros = FALSE, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays a histograms of a single asset

    # Arguments:
    #   x - a timeSeries object or any other rectangular object
    #       which can be transformed by the function as. matrix
    #       into a numeric matrix.

    # Example:
    #   x = as.timeSeries(data(LPP2005REC)) 
    #   par(mfrow = c(3,3)); assetsHistPlot(x); par(mfrow = c(1,1))
    
    # FUNCTION:

    # Settings:
    n = ncol(x)
    if (length(col) == 1) col = rep(col, times = n)

    # Plot:
    for (i in 1:n) {
        X = x[, i]
        if (skipZeros) X = X[series(X) != 0]
        histPlot(X, ylab = "Cumulated Returns", col = col[i], ...)
    }

    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


assetsLogDensityPlot =
    function(x, estimator = c("hubers", "sample", "both"),
    labels = TRUE, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays a pdf plot on logarithmic scale

    # Arguments:
    #   x - an uni- or multivariate return series of class 'timeSeries'
    #       or any other object which can be transformed by the function
    #       'as.timeSeries()' into an object of class 'timeSeries'.
    #   estimator - the type of estimator to fit the mean and variance
    #       of the density.
    #   doplot - a logical flag, by default TRUE. Should a plot be
    #       displayed?
    #   labels - a logical flag, by default TRUE. Should a default main
    #       title and labels addet to the plot?
    #   ... -

    # Details:
    #   Returns a pdf plot on a lin-log scale in comparison to a Gaussian
    #   density plot Two type of fits are available: a normal density with
    #   fitted sample mean and sample standard deviation, or a normal
    #   density with Hubers robust mean and standard deviation corfrected
    #   by the bandwidth of the Kernel estimator.

    # Example:
    #   x = as.timeSeries(data(LPP2005REC)) 
    #   par(mfrow=c(3,3)); assetsLogDensityPlot(x, "hubers"); par(mfrow=c(1,1))
    #   par(mfrow=c(3,3)); assetsLogDensityPlot(x, "sample"); par(mfrow=c(1,1))
    #   par(mfrow=c(3,3)); assetsLogDensityPlot(x, "both"); par(mfrow=c(1,1))
    
    # FUNCTION:

    # Settings:
    if (!is.timeSeries(x)) x = as.timeSeries(x)
    Units = colnames(x)
    doplot = TRUE

    # Select Type:
    estimator = match.arg(estimator)

    # Labels:
    if (labels) {
        main = "log PDF"
        xlab = "x"
        ylab = "log PDF"
    } else {
        main = xlab = ylab = ""
    }

    X = x

    for (i in 1:ncol(x)) {

        # Transform Data:
        x = as.vector(X[, i])
        if (labels) main = Units[i]

        # Kernel and Histogram Estimators:
        Density = density(x)
        Histogram = hist(x, breaks = "FD", plot = FALSE)
        result = list(density = Density, hist = Histogram)

        # Plot:
        if (doplot) {
            # Plot Frame:
            plot(Histogram$mids, log(Histogram$density), type = "n",
                lwd = 5, main = Units[i], xlab = xlab, ylab = ylab,
                xlim = range(Density$x), ylim = log(range(Density$y)),
                col = "red", ...)

            # Plot Density:
            points(Density$x, log(Density$y), pch = 19, col = "darkgrey",
                cex = 0.7)

            # Sample Line Fit:
            s = seq(min(Density$x), max(Density$x), length = 1001)
            if (estimator == "sample" || estimator == "both") {
                lines(s, log(dnorm(s, mean(x), sd(x))), col = "red", lwd = 2)
            }

            # Robust Huber Line Fit:
            if (estimator == "hubers" || estimator == "both") {
                h = MASS::hubers(x)
                logDensity = log(dnorm(s,
                    mean = h[[1]],
                    sd = sqrt(h[[2]]^2+Density$bw^2)))
                minLogDensity = log(min(Density$y))
                lines(
                    x = s[logDensity > minLogDensity],
                    y = logDensity[logDensity > minLogDensity],
                    col = "orange", lwd = 2)
            }

            # Plot Histogram:
            points(Histogram$mids, log(Histogram$density), pch = 19,
                col = "steelblue", ...)

            # Grid:
            if (labels) grid()
        }
    }

    # Return Value:
    invisible(result)
}


################################################################################

