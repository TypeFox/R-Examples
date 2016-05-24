
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
#  tFit                 Fits parameters of a Student-t density
################################################################################


tFit <-
function(x, df = 4, doplot = TRUE, span = "auto", trace = FALSE,
    title = NULL, description = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Return Maximum log-likelihood estimated

    # Note:
    #   Function Calls: nlminb(), density()

    # Example:
    #   tFit(rt(1000, df=4))

    # FUNCTION:

    # Transform:
    x.orig = x
    x = as.vector(x)

    # Settings:
    CALL = match.call()

    # Log-likelihood Function:
    obj = function(x, y = x, trace) {
        # Prevent from negative df's
        if (x[1] <= 0) x[1] = getRmetricsOptions(".x.save")
        f = -sum(log(dt(y, x[1])))
        # Print Iteration Path:
        if (trace) {
            cat("\n Objective Function Value:  ", -f)
            cat("\n Students df Estimate:      ", x[1], "\n")
        }
        setRmetricsOptions(.x.save = x[1])
        f
    }

    # Minimization:
    r = nlm(f = obj, p = c(df), y = x, trace = trace)

    # Optional Plot:
    if (doplot) {
        if (span == "auto") {
            df = r$estimate[1]
            span.min = qt(0.001, df)
            span.max = qt(0.999, df)
            span = seq(span.min, span.max, length = 100)
        }
        par(err = -1)
        z = density(x, n = 100, ...)
        x = z$x[z$y > 0]
        y = z$y[z$y > 0]
        y.points = dt(span, df = r$estimate[1])
        ylim = log(c(min(y.points), max(y.points)))
        plot(x, log(y), xlim = c(span[1], span[length(span)]),
            ylim = ylim, type = "p", xlab = "x", ylab = "log f(x)", ...)
        title("STUDENT-T: Parameter Estimation")
        lines(x = span, y = log(y.points), col = "steelblue")
        if (exists("grid")) grid()
    }

    # Add Title and Description:
    if (is.null(title)) title = "Student-t Parameter Estimation"
    if (is.null(description)) description = description()

    # Fit:
    fit = list(estimate = c(df = r$estimate), minimum = -r$minimum,
        code = r$code, gradient = r$gradient)

    # Return Value:
    new("fDISTFIT",
        call = as.call(CALL),
        model = "Student-t Distribution",
        data = as.data.frame(x.orig),
        fit = fit,
        title = as.character(title),
        description = description() )
}


################################################################################

