
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
# FUNCTION:                DESCRIPTION:
#  qqnormPlot               Returns a tailored normal quantile-quantile plot
#  qqnigPlot                Returns a tailored NIG quantile-quantile plot
#  qqghtPlot                Returns a tailored GHT quantile-quantile plot
#  qqgldPlot                Returns a tailored LD quantile-quantile plot
################################################################################


qqnormPlot <-
function(x, labels = TRUE, col = "steelblue", pch = 19,
    title = TRUE, mtext = TRUE, grid = FALSE, rug = TRUE, scale = TRUE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Example of a Normal quantile plot of data x to provide a visual
    #   assessment of its conformity with a normal (data is standardised
    #   first).

    # Arguments:
    #   x - an univariate return series of class 'timeSeries'
    #       or any other object which can be transformed by the function
    #       'as.timeSeries()' into an object of class 'timeSeries'.

    # Details:
    #   The ordered data values are posterior point estimates of the
    #   underlying quantile function. So, if you plot the ordered data
    #   values (y-axis) against the exact theoretical quantiles (x-axis),
    #   you get a scatter that should be close to a straight line if the
    #   data look like a random sample from the theoretical distribution.
    #   This function chooses the normal as the theory, to provide a
    #   graphical/visual assessment of how normal the data appear.
    #   To help with assessing the relevance of sampling variability on
    #   just "how close" to the normal the data appears, we add (very)
    #   approximate posterior 95% intervals for the uncertain quantile
    #   function at each point (Based on approximate theory) .

    # Author:
    #   Based on code written by Mike West, mw@stat.duke.edu

    # Note:
    #   Source from
    #   http://www.stat.duke.edu/courses/Fall99/sta290/Notes/

    # Example:
    #   x = rnorm(100); qqnormPlot(x); qqnormPlot(x, labels = FALSE)

    # FUNCTION:

    # Settings:
    if (!is.timeSeries(x)) x = as.timeSeries(x)
    Units = x@units
    x = as.vector(x)
    n = length(x)

    # Fit:
    p = (1:n)/(n+1)
    if (scale) x = (x-mean(x))/sqrt(var(x))
    par = c(mean = mean(x), var = var(x))

    # Quantiles:
    x = sort(x)
    p = ppoints(x)
    if (scale) z = qnorm(p) else z = qnorm(p, mean(x), sd(x))

    # Plot:
    if (labels) {
        xlab = "Normal Quantiles"
        ylab = paste(Units, "Ordered Data")
        plot(z, x, xlab = xlab, ylab = ylab,
            col = col, pch = 19, ...)
    } else {
        plot(z, x, ...)
    }

    # Title:
    if(title) {
        title(main = "NORM QQ PLOT")
    }

    # Margin Text:
    if (mtext) {
        Text = "Confidence Intervals: 95%"
        mtext(Text, side = 4, adj = 0, col = "darkgrey", cex = 0.7)
    }

    # Grid:
    if (grid) {
        grid()
    }

    # Add Diagonal Line:
    abline(0, 1, col = "grey")

    # Add Rugs:
    if(rug) {
        rug(z, ticksize = 0.01, side = 1, quiet = TRUE)
        rug(x, ticksize = 0.01, side = 2, quiet = TRUE)
    }

    # 95% Confidence Intervals:
    s = 1.96*sqrt(p*(1-p)/n)
    pl = p-s
    i = pl<1 & pl>0
    lower = quantile(x, probs = pl[i])
    lines(z[i], lower, col = "brown")
    pl = p+s
    i = pl < 1 & pl > 0
    upper = quantile(x, probs = pl[i])
    lines(z[i], upper, col = "brown")
    abline(h = mean(x), col = "grey")
    abline(v = mean(x), col = "grey")

    # Result:
    ans = list(x = z, y = x)
    attr(ans, "control")<-par

    # Return Value:
    invisible(ans)
}


# ------------------------------------------------------------------------------


qqnigPlot <-
function(x, labels = TRUE, col = "steelblue", pch = 19,
    title = TRUE, mtext = TRUE, grid = FALSE, rug = TRUE, scale = TRUE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays a NIG quantile-quantile Plot

    # Arguments:
    #   x - an univariate return series of class 'timeSeries'
    #       or any other object which can be transformed by the function
    #       'as.timeSeries()' into an object of class 'timeSeries'.

    # Example:
    #   qqnigPlot(rnig(100))

    # FUNCTION:

    # Settings:
    if (!is.timeSeries(x)) x = as.timeSeries(x)
    stopifnot(isUnivariate(x))
    Units = x@units
    x = as.vector(x)
    n = length(x)

    ## YC: no scaling
    ## FIXME: should take care of too small time series

    # Fit:
    fit = nigFit(x, doplot = FALSE, trace = FALSE)
    par = fit@fit$estimate
    names(par) = c("alpha", "beta", "delta", "mu")

    # Quantiles:
    x = sort(x)
    p = ppoints(x)

    z = qnig(p, par[1], par[2], par[3], par[4])

    # Plot:
    if (labels) {
        xlab = "Theoretical Quantiles"
        ylab = "Sample Quantiles"
        plot(z, x, xlab = xlab, ylab = ylab, col = col, pch = pch, ...)
    } else {
        plot(z, x, ...)
    }

    # Title:
    if (title) {
        title(main = "NIG QQ Plot")
    }

    # Margin Text:
    rpar = signif(par, 3)
    text = paste(
        "alpha =", rpar[1],
        "| beta =", rpar[2],
        "| delta =", rpar[3],
        "| mu =", rpar[4])
    mtext(text, side = 4, adj = 0, col = "grey", cex = 0.7)

    # Grid:
    if (grid) {
        grid()
    }

    # Add Fit:
    abline(lsfit(z, x))

    # Add Rugs:
    if(rug) {
        rug(z, ticksize = 0.01, side = 3, quiet = TRUE)
        rug(x, ticksize = 0.01, side = 4, quiet = TRUE)
    }

    # Result:
    ans = list(x = z, y = x)
    attr(ans, "control") <- par

    # Return Value:
    invisible(ans)
}


# ------------------------------------------------------------------------------


qqghtPlot <-
function(x, labels = TRUE, col = "steelblue", pch = 19,
    title = TRUE, mtext = TRUE, grid = FALSE, rug = TRUE, scale = TRUE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays a GHT quantile-quantile Plot

    # Arguments:
    #   x - an univariate return series of class 'timeSeries'
    #       or any other object which can be transformed by the function
    #       'as.timeSeries()' into an object of class 'timeSeries'.

    # Example:
    #   qqnigPlot(rgh(100))

    # FUNCTION:

    # Settings:
    if (!is.timeSeries(x)) x = as.timeSeries(x)
    stopifnot(isUnivariate(x))
    Units = x@units
    x = as.vector(x)
    n = length(x)

    # Fit:
    fit = ghtFit(x, doplot = FALSE, trace = FALSE)
    par = fit@fit$estimate
    names(par) = c("beta", "delta", "mu", "nu")

    # Plot:
    x <- sort(x)
    p <- ppoints(x)
    z <- qght(p, par[1], par[2], par[3], par[4])


    if (labels) {
        plot(z, x, col = col, ann = FALSE, ...)
    } else {
        plot(z, x, ...)
    }

    # Add Grid:
    if (grid) {
        grid()
    }

    # Add title:
    if (title) {
        title(
            main = "GHT QQ Plot",
            xlab = "Theoretical Quantiles",
            ylab = "Sample Quantiles")
    }


    # Add Fit:
    abline(lsfit(z, x))

    # Add Rugs:
    if(rug) {
        rug(z, ticksize = 0.01, side = 3, quiet = TRUE)
        rug(x, ticksize = 0.01, side = 4, quiet = TRUE)
    }

    # Result:
    ans = list(x = z, y = x)
    attr(ans, "control")<-par

    # Return Value:
    invisible(ans)
}


################################################################################


qqgldPlot <-
function(x, labels = TRUE, col = "steelblue", pch = 19,
    title = TRUE, mtext = TRUE, grid = FALSE, rug = TRUE, scale = TRUE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays a Generalized lambda Distribution quantile-quantile Plot

    # Arguments:
    #   x - an univariate return series of class 'timeSeries'
    #       or any other object which can be transformed by the function
    #       'as.timeSeries()' into an object of class 'timeSeries'.

    # Example:
    #   qqgldPlot(rgld(100))

    # FUNCTION:

    # Settings:
    if (!is.timeSeries(x)) x = as.timeSeries(x)
    stopifnot(isUnivariate(x))
    Units = x@units
    x = as.vector(x)
    n = length(x)

    ## YC: no scaling
    ## FIXME: should take care of too small time series

    # Fit:
    fit = gldFit(x, doplot = FALSE, trace = FALSE)
    par = fit@fit$estimate
    names(par) = c("lambda1", "lambda2", "lambda3", "lambda4")

    # Quantiles:
    x = sort(x)
    p = ppoints(x)

    z = qgld(p, par[1], par[2], par[3], par[4])

    # Plot:
    if (labels) {
        xlab = "Theoretical Quantiles"
        ylab = "Sample Quantiles"
        plot(z, x, xlab = xlab, ylab = ylab, col = col, pch = pch, ...)
    } else {
        plot(z, x, ...)
    }

    # Title:
    if (title) {
        title(main = "NIG QQ Plot")
    }

    # Margin Text:
    rpar = signif(par, 3)
    text = paste(
        "lambda1 =", rpar[1],
        "| lambda2 =", rpar[2],
        "| lambda3 =", rpar[3],
        "| lambda4 =", rpar[4])
    mtext(text, side = 4, adj = 0, col = "grey", cex = 0.7)

    # Grid:
    if (grid) {
        grid()
    }

    # Add Fit:
    abline(lsfit(z, x))

    # Add Rugs:
    if(rug) {
        rug(z, ticksize = 0.01, side = 3, quiet = TRUE)
        rug(x, ticksize = 0.01, side = 4, quiet = TRUE)
    }

    # Result:
    ans = list(x = z, y = x)
    attr(ans, "control") <- par

    # Return Value:
    invisible(ans)
}


################################################################################