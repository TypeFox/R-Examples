
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
# FUNCTION:            DESCRIPTION:
#  .normFit             Fits parameters of a Normal density
################################################################################


# normFit is now in fBasics


# ------------------------------------------------------------------------------


.normFit <-
function(x, doplot = TRUE, span = "auto", title = NULL,
    description = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Return Maximum log-likelihood estimated
    #   Paramters for Normal Distribution

    # Notes:
    #   Function Calls: nlminb(), density()
    #   The function normFit() can be found in the Rmetrics
    #       chapter GarchDistributions.

    # FUNCTION:

    # Transform:
    x.orig = x
    x = as.vector(x)

    # Settings:
    CALL = match.call()

    # MLE:
    N = length(x)
    mean = sum(x)/N
    sd = sqrt(sum((x-mean)^2)/N)

    # Optional Plot:
    if (doplot) {
        if (span == "auto") {
            span.min = qnorm(0.001, mean, sd)
            span.max = qnorm(0.999, mean, sd)
            span = seq(span.min, span.max, length = 100)
        }
        par(err = -1)
        z = density(x, n = 100, ...)
        x = z$x[z$y > 0]
        y = z$y[z$y > 0]
        y.points = dnorm(span, mean, sd)
        ylim = log(c(min(y.points), max(y.points)))
        plot(x, log(y), xlim = c(span[1], span[length(span)]),
            ylim = ylim, type = "p", xlab = "x", ylab = "log f(x)", ...)
        title("NORMAL: Parameter Estimation")
        lines(x = span, y = log(y.points), col = "steelblue")
        if (exists("grid")) grid()
    }

    # Add Title and Description:
    if (is.null(title)) title = "Normal Parameter Estimation"
    if (is.null(description)) description = description()

    # Fit:
    fit = list(
        estimate = c(mean = mean, sd = sd),
        minimum = sum(log(dnorm(x, mean, sd))),
        code = NA)

    # Return Value:
    new("fDISTFIT",
        call = as.call(CALL),
        model = "Normal Distribution",
        data = as.data.frame(x.orig),
        fit = fit,
        title = as.character(title),
        description = description() )
}


# ------------------------------------------------------------------------------


nFit <- .normFit


################################################################################

