
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


################################################################################
# FUNCTION:             DESCRIPTION:
#  scalinglawPlot        Evaluates and displays scaling law behavior
################################################################################


scalinglawPlot <- 
function(x, span = ceiling(log(length(x)/252)/log(2)), doplot = TRUE,
    labels = TRUE, trace = TRUE, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Investigates the scaling law.
    #   The input "x" requires log-returns.

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

    # Labels:
    if (labels) {
        main = "Scaling Law Plot"
        xlab = "log-time"
        ylab = "log-volatility"
    } else {
        main = xlab = ylab = ""
    }

    X = x
    DIM = dim(X)[2]
    Intercept = Exponent = InverseExponent = NULL
    for (i in 1:DIM) {

        # Get Data:
        x = as.vector(as.matrix(X)[, i])
        if (labels) main = Units[i]


        # Settings:
        logtimesteps = span
        xmean = mean(x)

        # x have to be logarithmic returns
        y = (x-xmean)
        logprices = cumsum(y)

        # Scaling Power Low:
        scale = function(nx, logprices) {
            sum(abs(diff(logprices, lag = (2^nx))))}
        nx = 0:logtimesteps; x = nx*log(2)
        y = log(apply(matrix(nx), 1, scale, logprices))
        # fit = lsfit(x, y)$coefficients

        # Runs in both environments, R and SPlus:
        fit = lsfit(x, y)

        # Robust Fit:
        # fit = l1fit(x, y)

        # Fit Result:
        Fit = unlist(fit)[1:2]
        alpha = 1.0/Fit[2]
        if (doplot) {
            plot(x, y, main = main, xlab = xlab, ylab = ylab, ...)
            abline(Fit[1], Fit[2])
            abline(Fit[1], 0.5, col = "steelblue")
        }
        if (labels) grid()

        # Trace:
        if (trace) {
            cat ("\nScaling Law:         ", Units[i])
            cat ("\n  Plot Intercept     ", fit$coefficients[1])
            cat ("\n  Plot Slope         ", fit$coefficients[2])
            cat ("\n  Plot Inverse Slope ", 1/fit$coefficients[2])
            cat ("\n\n")
        }
        Intercept = c(Intercept, fit$coefficients[1])
        Exponent = c(Exponent, fit$coefficients[2])
        InverseExponent = c(InverseExponent, 1/fit$coefficients[2])

    }

    names(Intercept) = Units
    names(Exponent) = Units
    names(InverseExponent) = Units
    result = list(
        Intercept = Intercept,
        Exponent = Exponent,
        InverseExponent = InverseExponent)

    # Return Value:
    invisible(result)
}


################################################################################

