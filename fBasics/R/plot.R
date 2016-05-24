
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
# FUNCTION:                 INTERNAL USED PLOT FUNCTIONS:
#  .residualsPlot            Returns a residual series plot
#  .acfPlot                  Returns a autocorrelation function plot
#  .pacfPlot                 Returns a partial ACF plot
#  .mrlPlot                  Returns a mean residual life plot
# FUNCTION:                 INTERNAL USED BIVARIATE PLOT FUNCTIONS:
#  .responsesPlot            Returns a response series plot
#  .firePlot                 Returns a fitted values vs.residuals plot
# FUNCTION:                 INTERNAL THREE-DIMENSIONAL PLOT UTILITIES:
#  .circlesPlot              Returns a circles plot indexing 3rd variable
#  .perspPlot                Returns a perspective plot in 2 dimensions
#  .contourPlot              Returns a contour plot in 2 dimensions
#  .histStack                Returns a stacked histogram plot
################################################################################


.residualsPlot <- 
function(x, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns time series graph of residuals

    # Arguments:
    #   x - an univariate time series of residuals

    # FUNCTION:

    # Get Data:
    x = as.vector(x)

    # Plot:
    plot(x, type = "l", ylab = "Residuals",
        main = "Residual Series", col = "steelblue", ...)
    rug(x, ticksize = 0.01, side = 4)
    grid()
    abline(h = 0, col = "grey")

    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


.acfPlot <- 
function(x, ...)
{
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION:

    # Convert Type:
    x = as.vector(x)

    # ACF:
    acf(x, ...)

    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


.pacfPlot <- 
function(x, ...)
{
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION:

    # Convert Type:
    x = as.vector(x)

    # ACF:
    pacf(x, ...)

    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


.mrlPlot <- 
function(x, ci = 0.95, umin = mean(x), umax = max(x), nint = 100,
    doplot = TRUE, plottype = c("autoscale", ""), labels = TRUE, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Create a mean residual life plot with
    #   confidence intervals.

    # References:
    #   A function originally written by S. Coles

    # FUNCTION:

    # Convert Type:
    x = as.vector(x)

    # Settings:
    plottype = plottype[1]

    # Convert x to a vector, if the input is a data.frame.
    if (is.data.frame(x)) x = x[,1]
    sx = xu = xl = rep(NA, nint)
    u = seq(umin, umax, length = nint)
    for (i in 1:nint) {
        x = x[x >= u[i]]
        sx[i] = mean(x - u[i])
        sdev = sqrt(var(x))
        n = length(x)
        xu[i] = sx[i] + (qnorm((1 + ci)/2) * sdev) / sqrt(n)
        xl[i] = sx[i] - (qnorm((1 + ci)/2) * sdev) / sqrt(n)
    }

    # Plot:
    if (doplot) {
        if (labels) {
            xlab = "Threshold: u"
            ylab = "Mean Excess: e"
            main = "Mean Residual Live Plot"
        } else {
            main = xlab = ylab = ""
        }
        if (plottype == "autoscale") {
            ylim = c(min(xl[!is.na(xl)]), max(xu[!is.na(xu)]))
            plot(u, sx, type = "o", pch = 19, col = "steelblue",
                xlab = xlab, ylab = ylab, ylim = ylim, main = main, ...)
        } else {
            plot(u[!is.na(xl)], sx[!is.na(xl)], type = "o",
                pch = 19, col = "steelblue",
                xlab = xlab, ylab = ylab, main = main, ...)
        }
        lines(u[!is.na(xl)], xl[!is.na(xl)], col = "brown")
        lines(u[!is.na(xu)], xu[!is.na(xu)], col = "brown")
        if (labels) {
            grid()
            text = paste("ci =", as.character(round(ci, 3)))
            mtext(text, side = 4, adj = 0, cex = 0.7)
        }
    }

    # Result
    result = data.frame(threshold = u, mrl = sx)

    # Return Value:
    if (doplot) return(invisible(result)) else return(result)
}


################################################################################


.responsesPlot <- 
function(x, y = NULL, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns time series graph of responses and fitted values

    # Arguments:
    #   x - an univariate time series of responses
    #   y - an univariate time series of fitted values

    # FUNCTION:

    # Get Data:
    x = as.vector(x)
    y = as.vector(y)

    # Responses Plot:
    plot(x, type = "l", ylab = "Responses",
        main = "Responses & Fitted Values", col = "steelblue", ...)
    rug(x, ticksize = 0.01, side = 4)
    grid()
    abline(h = 0, col = "grey")

    # Add fitted values:
    if (!is.null(y)) points(y, pch = 19, col = "red")

    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


.firePlot <- 
function(x, y, method = c("scatter", "hist"), ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns fitted values vs. residuals plots

    # Arguments:
    #   x - univariate time series (residuals)
    #   y - univariate time series (fitted)

    # FUNCTION:

    # Check Arguments:
    method = match.arg(method)

    # Get Data:
    x = as.vector(x)
    y = as.vector(y)


    if (method == "scatter") {

        # Scatter Plot:
        plot(x, y,
            xlab = "Fitted Values", ylab = "Residuals",
            main = "Residuals vs Fitted",
            pch = 19, col = "steelblue")
        panel.smooth(x, y)
        abline(h = 0, lty = 3, col = "grey")
        rug(x, ticksize = 0.01, side = 3)
        rug(y, ticksize = 0.01, side = 4)

    } else if (method == "hist") {

        # Histogram Plot:

        # Save default, for resetting ...
        def.par = par(no.readonly = TRUE)

        # Layout:
        nf = layout(matrix(c(2, 0, 1, 3), 2, 2, byrow = TRUE), c(3, 1),
            c(1, 3), TRUE)

        # Scatterplot:
        par(mar = c(3 ,3, 1, 1))
        plot(x, y, xlim = range(x), ylim = range(y),
            xlab = "", ylab = "", pch = 19, col = "steelblue")
        panel.smooth(x, y)
        abline(h = 0, lty = 3, col = "grey")
        rug(x, side = 3)
        rug(y, side = 4)

        # Histogram:
        xhist = hist(x, nclass = 15, plot = FALSE)
        yhist = hist(y, nclass = 15, plot = FALSE)
        top = max(c(xhist$counts, yhist$counts))

        # X-Side:
        par(mar = c(0, 3, 1, 1))
        Main = "\n                            Fitted"
        barplot(xhist$counts, axes = FALSE, ylim = c(0, top),
            space = 0, col = "steelblue", border = "white",
            main = Main)
        abline(h = 0, lwd = 2, col = "grey")

        # Y-Side:
        par(mar = c(3, 0, 1, 1))
        barplot(yhist$counts, axes = FALSE, xlim = c(0, top),
            space = 0, col = "steelblue", , border = "white",
            horiz = TRUE, main = "Residuals")
        abline(v = 0, lwd = 2, col = "grey")

        # Reset:
        par(def.par)

    }

    # Return Value:
    invisible()
}


################################################################################


.circlesPlot <- 
function(x, y = NULL, z = NULL, scale = 1, points = TRUE,
    labels = TRUE, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Creates a scatterplot with circle size as third variable

    # Example:
    #   circlesPlot(x=rnorm(50), y=rnorm(50), z=rnorm(50))
    #   circlesPlot(x=rnorm(50), y=rnorm(50), z=rnorm(50), labels= FALSE)

    # FUNCTION:

    # Transfor Input:
    if (is.list(x)) x = matrix(unlist(x), ncol = 3)
    if (is.data.frame(x)) x = as.matrix.data.frame(x)
    if (is.matrix(x)) {
        z = x[, 3]
        y = x[, 2]
        x = x[, 1]
    }
    nX = length(x)
    nY = length(y)
    # nZ = length(z)
    stopifnot(nX == nY)
    # stopifnot(nX == nZ || nX*nY == nZ)

    # Create Circle Plot:
    if (labels) {
        plot(x, y, type = "n")
    } else {
        plot(x, y, xlab = "", ylab = "", type = "n")
    }
    symbols(x, y, add = TRUE, circles = abs(z)^scale, inches = 0.25,
        fg = "black", bg = "steelblue", ...)
    X = x[z < 0]
    Y = y[z < 0]
    Z = z[z < 0]
    symbols(X, Y, add = TRUE, circles = abs(Z)^scale, inches = 0.25,
        fg = "black", bg = "orange", ...)
    if (points) points(x, y, pch = 19)
    grid()

    # Return Value:
    invisible(NULL)
}


# ------------------------------------------------------------------------------


.perspPlot <- 
function(x, y, z, theta = -40, phi = 30, col = "steelblue", ps = 9, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns a perspecvtive plot

    # Notes:
    #   A synonyme call for function 'persp'

    # FUNCTION:

    # Perspective Plot:
    if (class(version) == "Sversion") {
        # we assume SPlus:
        ans = persp(x = x, y = y, z = z, ...)
    } else {
        # R:
        par(ps = ps)
        if (!exists("ticktype")) ticktype = "detailed"
        if (!exists("expand")) expand = 0.6
        if (!exists("r")) r = 500
        ans = persp(x = x, y = y, z = z, theta = theta, phi = phi,
            col = col, ticktype = ticktype, expand = expand, ...)
    }

    # Return Value:
    invisible(ans)
}


# ------------------------------------------------------------------------------


.contourPlot <- 
function(x, y, z, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns a contour plot

    # Notes:
    #   A synonyme call for function 'contour'

    # FUNCTION:

    # Contour Plot:
    if (class(version) == "Sversion") {
        # we assume SPlus:
        ans = contour(x = x, y = y, z = z, ...)
    } else {
        # R:
        ans = contour(x = x, y = y, z = z, ...)
    }

    # Return Value:
    invisible(ans)
}


# ------------------------------------------------------------------------------


.histStack <- 
function(x, y = NULL, space = 0, ylab = "frequency", ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns a stacked histogram Plot

    # Example:
    #   .histStack(rnorm(1000, -1), rnorm(1000, 1))

    # FUNCTION:

    # Compute Histograms:
    breaks = hist(c(x, y))$breaks
    bars = rbind(
        hist(x, breaks = breaks, plot = FALSE)$counts,
        hist(y, breaks = breaks, plot = FALSE)$counts)

    # Plot:
    barplot(bars, space = space, ylab = ylab, ...)
    at = seq(along = breaks) - 1
    modulo = ceiling(length(at)/10)
    sel = (at%%modulo == 0)
    axis(side = 1, at = at[sel], labels = paste(breaks)[sel])

    # Return Value:
    invisible()
}


################################################################################

