
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
# FUNCTION:            GENERALIZED DISTRIBUTION:
#  ghFit                Fits parameters of a generalized hyperbolic DF
################################################################################


ghFit <-
function(x, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = -1/2,
    scale = TRUE, doplot = TRUE, span = "auto", trace = TRUE,
    title = NULL, description = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits parameters of a generalized hyperbolic density

    # FUNCTION:

    # Transform:
    x.orig = x
    x = as.vector(x)
    if (scale) {
        SD = sd(x)
        x = x / SD
    }

    # Settings:
    CALL = match.call()

    # Log-likelihood Function:
    obj = function(x, y = x, trace){
        if (NA %in% x) return(1e99)
        if (abs(x[2]) >= x[1]) return(1e99)
        f = -sum(dgh(y, x[1], x[2], x[3], x[4], x[5], log = TRUE))
        # Print Iteration Path:
        if (trace) {
            cat("\n Objective Function Value:  ", -f)
            cat("\n Parameter Estimates:       ", x, "\n")
        }
        f
    }

    # Minimization:
    r = # Variable Transformation and Minimization:
    eps = 1e-10
    BIG = 1000
    f = obj(x = c(alpha, beta, delta, mu, lambda), y = x, trace = FALSE)
    r = nlminb(start = c(alpha, beta, delta, mu, lambda), objective = obj,
        lower = c(eps, -BIG, eps, -BIG, -BIG), upper = BIG, y = x,
        trace = trace)
    names(r$par) <- c("alpha", "beta", "delta", "mu", "lambda")

    # Result:
    if (scale) {
        r$par = r$par / c(SD, SD, 1/SD, 1/SD, 1)
        r$objective = obj(r$par, y = as.vector(x.orig), trace = trace)
    }

    # Optional Plot:
    if (doplot) {
        x = as.vector(x.orig)
        if (span == "auto") span = seq(min(x), max(x), length = 51)
        z = density(x, n = 100, ...)
        x = z$x[z$y > 0]
        y = z$y[z$y > 0]
        y.points = dnig(span, r$par[1], r$par[2], r$par[3], r$par[4])
        ylim = log(c(min(y.points), max(y.points)))
        plot(x, log(y), xlim = c(span[1], span[length(span)]),
            ylim = ylim, type = "p", xlab = "x", ylab = "log f(x)", ...)
        title("GH Parameter Estimation")
        lines(x = span, y = log(y.points), col = "steelblue")
    }

    # Add Title and Description:
    if (is.null(title)) title = "Generalized Hyperbolic Parameter Estimation"
    if (is.null(description)) description = description()

    # Fit:
    fit = list(estimate = r$par, minimum = -r$objective, code = r$convergence)

    # Return Value:
    new("fDISTFIT",
        call = as.call(CALL),
        model = "Generalized Hyperbolic Distribution",
        data = as.data.frame(x.orig),
        fit = fit,
        title = as.character(title),
        description = description() )
}


################################################################################

