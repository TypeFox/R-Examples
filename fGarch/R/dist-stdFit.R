
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
# FUNCTION:              PARAMETER ESTIMATION:
#  stdFit                 Fit the parameters for a Sudent-t distribution
################################################################################


.stdFit <- 
function(x, mean = 0, sd = 1, xi = 1.5, 
    scale = NA, doplot = TRUE, add = FALSE, span = "auto", trace = TRUE,
    title = NULL, description = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits parameters of a Student-t using maximum log-likelihood  

    # Example:
    #   set.seed(4711); x = rstd(500); .stdFit(x)@fit$estimate 
    
    # FUNCTION:

    # Settings:
    dist = dstd
    model = "STD Parameter Estimation"
    scale = "not used"
    x = x.orig = as.vector(x)

    # Parameter Estimation:
    obj = function(x, y = x, trace) {
        f = try(-sum(log(dist(y, x[1], x[2], x[3]))), silent = TRUE)
        if (is.na(f) | class(f) == "try-error") return(1e9)
        # Print Iteration Path:
        if (trace) {
            cat("\n Objective Function Value:  ", -f)
            cat("\n Parameter Estimates:       ", x, "\n")
        }
        f }
    r = nlminb(
        start = c(mean = 0, sd = 1, nu = 5), 
        objective = obj,
        lower = c(-Inf, 0,   2), 
        upper = c( Inf, Inf, Inf), 
        y = x, 
        trace = trace)
    names(r$par) <- c("mean", "sd", "nu")

    # Add Title and Description:
    if (is.null(title)) title = model
    if (is.null(description)) description = description()

    # Result:
    fit = list(estimate = r$par, minimum = -r$objective, code = r$convergence)

    # Optional Plot:
    if (doplot) {
        x = as.vector(x.orig)
        if (span == "auto") span = seq(min(x), max(x), length = 501)
        z = density(x, n = 100, ...)
        x = z$x[z$y > 0]
        y = z$y[z$y > 0]
        y.points = dstd(span, r$par[1], r$par[2], r$par[3])
        ylim = log(c(min(y.points), max(y.points)))
        if (add) {
            lines(x = span, y = log(y.points), col = "steelblue")
        } else {
            plot(x, log(y), xlim = c(span[1], span[length(span)]),
                ylim = ylim, type = "p", xlab = "x", ylab = "log f(x)", ...)
            title(main = model)
            lines(x = span, y = log(y.points), col = "steelblue")
        }
    }

    # Return Value:
    new("fDISTFIT",
        call = match.call(),
        model = model,
        data = as.data.frame(x.orig),
        fit = fit,
        title = title,
        description = description() )
}


# ------------------------------------------------------------------------------


stdFit <-
function(x, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Fit the parameters for a skew Normal distribution

    # FUNCTION:

    # Start Value:
    start = c(mean = mean(x), sd = sqrt(var(x)), nu = 4)

    # Log-likelihood Function:
    loglik = function(x, y = x){
        f = -sum(log(dstd(y, x[1], x[2], x[3])))
        f }

    # Minimization:
    fit = nlminb(
        start = start, 
        objective = loglik,
        lower = c(-Inf, 0, 2), 
        upper = c(Inf, Inf, Inf), 
        y = x, ...)

    # Add Names to $par
    names(fit$par) = c("mean", "sd", "nu")

    # Return Value:
    fit
}


################################################################################

