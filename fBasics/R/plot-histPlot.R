
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


################################################################################
# FUNCTION:                TAILORED DENSITY FUNCTIONS:
#  histPlot                 Returns a tailored histogram plot
#  densityPlot              Returns a tailored kernel density estimate plot
#  logDensityPlot           Returns a tailored log kernel density estimate plot
#  .plot.histogram          Replaces here the function plot.histogram
################################################################################


histPlot <-
function(x, labels = TRUE, col = "steelblue", fit = TRUE,
    title = TRUE, grid = TRUE, rug = TRUE, skip = FALSE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns a probability histogram plot for each column of a
    #   timeSeries object

    # Arguments:
    #   x - an uni- or multivariate return series of class 'timeSeries'
    #       or any other object which can be transformed by the function
    #       'as.timeSeries()' into an object of class 'timeSeries'.

    # FUNCTION:

    # timeSeries:
    stopifnot(is.timeSeries(x))
    N = NCOL(x)
    Units = colnames(x)
    if (length(col) == 1) col = rep(col, times = N)

    # Histogram Plots:
    for (i in 1:N) {

        # Histogram:
        X = as.vector(x[, i])
        if (skip) X = X[X != 0]

        # Plot:
        if (labels) {
            H = hist(x = X, , breaks = "FD", plot = FALSE)
            .plot.histogram(H, col = col[i], freq = FALSE, ...)
            box()
        } else {
            H = hist(x = X, plot = FALSE, ...)
            .plot.histogram(H, col = col[i], freq = FALSE, ...)
        }

        # Add Title:
        if(title) {
            title(main = paste(Units[i], "Histogram"),
                xlab = "Value", ylab = "Probability")
        }

        # Add Fit:
        if (fit) {
            mean = mean(X)
            sd = sd(X)
            xlim = range(H$breaks)
            s = seq(xlim[1], xlim[2], length = 201)
            lines(s, dnorm(s, mean, sd), lwd = 2, col = "brown")
        }

        # Add Mean:
        if (labels) {
            abline(v = mean, lwd = 2, col = "orange")
            Text = paste("Mean:", signif(mean, 3))
            mtext(Text, side = 4, adj = 0, col = "darkgrey", cex = 0.7)
        }

        # Add Grid:
        if (grid) grid()

        # Add Zero Line:
        if(labels) {
            abline(h = 0, col = "grey")
        }

        # Add Rug Plot:
        if(rug) {
            rug(X, ticksize = 0.01, quiet = TRUE)
        }
    }

    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


densityPlot <-
function(x, labels = TRUE, col = "steelblue", fit = TRUE, hist = TRUE,
    title = TRUE, grid = TRUE, rug = TRUE, skip = FALSE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns density plots for each column of a timeSeries object

    # FUNCTION:

    # timeSeries:
    stopifnot(is.timeSeries(x))
    N = NCOL(x)
    Units = colnames(x)
    if (length(col) == 1) col = rep(col, times = N)

    # Density Plots:
    for (i in 1:N) {

        # Density:
        X = as.vector(x[, i])
        if (skip) X = X[X != 0]

        # Underlay Histogram:
        if (hist) {
            H = hist(x = X, , breaks = "FD", plot = FALSE)
            plot(x = H$mids, y = H$density, type = "h", lwd = 2,
                main = "", xlab = "", ylab = "", col = "grey")
        }

        # Plot Density:
        D = density(X, ...)
        if (hist) {
            lines(D$x, D$y, lwd = 2, col = "brown")
        } else {
            plot(D, col = col[i], ann = FALSE, ...)
        }

        # Add Title:
        if (title) {
            title(main = Units[i], xlab = "Value", ylab = "Density")
        }

        # Add Fit:
        if (fit) {
            mean = mean(X)
            sd = sd(X)
            xlim = range(H$breaks)
            s = seq(xlim[1], xlim[2], length = 201)
            lines(s, dnorm(s, mean, sd), lwd = 2, col = "darkgrey")
        }

        # Add Mean:
        if (labels) {
            abline(v = mean, lwd = 2, col = "orange")
            Text = paste("Mean:", signif(mean, 3))
            mtext(Text, side = 4, adj = 0, col = "darkgrey", cex = 0.7)
        }

        # Add Grid:
        if (grid) {
            grid()
        }

        # Add Zero Line:
        if(labels) {
            abline(h = 0, col = "grey")
        }

        # Add Rug Plot:
        if(rug) {
            rug(X, ticksize = 0.01, quiet = TRUE)
        }
    }

    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


logDensityPlot <-
function(x, labels = TRUE, col = "steelblue", robust = TRUE,
    title = TRUE, grid = TRUE, rug = TRUE, skip = FALSE,  ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays a pdf plot on logarithmic scale

    # Details:
    #   Returns a pdf plot on a lin-log scale in comparison to a Gaussian
    #   density plot Two type of fits are available: a normal density with
    #   fitted sample mean and sample standard deviation, or a normal
    #   density with Hubers robust mean and standard deviation corfrected
    #   by the bandwidth of the Kernel estimator.

    # FUNCTION:

    # timeSeries:
    stopifnot(is.timeSeries(x))
    N = NCOL(x)
    Units = colnames(x)
    if (length(col) == 1) col = rep(col, times = N)

    # Log Density:
    for (i in 1:N) {

        # Transform Data:
        X = as.vector(x[, i])
        if (skip) X = X[X != 0]

        # Kernel and Histogram Estimators:
        Density = density(X, ...)
        Histogram = hist(X, breaks = "FD", plot = FALSE)

        # Plot Frame:
        plot(Histogram$mids, log(Histogram$density), type = "n",
            ann = FALSE,
            xlim = range(Density$x), ylim = log(range(Density$y)), ...)

        # Add Title:
        if(title) {
            title(main = paste(Units[i], "Log Density"),
                xlab = "Value", ylab = "Log Density")
        }

        # Add Kernel Density Estimated Points:
        points(Density$x, log(Density$y), pch = 19, cex = 0.7, col = "grey")

        # Sample Line Fit:
        s = seq(min(Density$x), max(Density$x), length = 1001)

        # Robust Fit:
        if (robust) {
            h = MASS::hubers(X)
            logDensity = log(dnorm(s,
                mean = h[[1]],
                sd = sqrt(h[[2]]^2+Density$bw^2)))
            minLogDensity = log(min(Density$y))
            lines(
                x = s[logDensity > minLogDensity],
                y = logDensity[logDensity > minLogDensity],
                col = "red", lwd = 2)

            # Standard Fit:
            lines(s, log(dnorm(s, mean(X), sd(X))),
                col = "orange", lwd = 2)
        }

        # Plot Histogram:
        points(Histogram$mids, log(Histogram$density),
            pch = 151, col = "steelblue")

        # Grid:
        if (grid) {
            grid()
        }

        # Add Rug Plot:
        if(rug) {
            rug(x, ticksize = 0.01, quiet = TRUE)
        }
    }

    # Return Value:
    invisible()
}


################################################################################


.plot.histogram <-
function(x, freq = equidist, density = NULL, angle = 45,
    col = NULL, border = "white", lty = NULL,
    main = NULL, sub = NULL, xlab = NULL, ylab = NULL,
    xlim = range(x$breaks), ylim = NULL,
    axes = TRUE, labels = FALSE, add = FALSE, ...)
{
    # This replacement of plot.histogram() suppresses title
    #   printing which would be otherwise not possible!

    equidist <- if(is.logical(x$equidist)) x$equidist
    else { h <- diff(x$breaks) ; diff(range(h)) < 1e-7 * mean(h) }
    if(freq && !equidist)
    warning("the AREAS in the plot are wrong -- rather use freq=FALSE")

    y <- if (freq) x$counts else { ## x$density -- would be enough, but
    ## for back compatibility
    y <- x$density; if(is.null(y)) x$intensities else y}
    nB <- length(x$breaks)
    if(is.null(y) || 0 == nB) stop("'x' is wrongly structured")

    if(!add) {
        if(is.null(ylim)) ylim <- range(y, 0)
        if (missing(ylab)) ylab <- if (!freq) "Density" else "Frequency"
        plot.new()
        plot.window(xlim, ylim, "") #-> ylim's default from 'y'
        # title(main = main, sub = sub, xlab = xlab, ylab = ylab, ...)
        if(axes) { axis(1, ...)
            axis(2, ...) } }
    rect(x$breaks[-nB], 0, x$breaks[-1], y,
        col = col, border = border,
        angle = angle, density = density, lty = lty)
    if((logl <- is.logical(labels) && labels) || is.character(labels))
    text(x$mids, y,
        labels = if(logl) { if(freq) x$counts else round(x$density,3)
        } else labels, adj = c(0.5, -0.5))
    invisible()
}


################################################################################

