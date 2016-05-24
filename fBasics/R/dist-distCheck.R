
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
# FUNCTION:                 DESCRIPTION:
#  distCheck                Checks consistency of distributions
################################################################################


distCheck <- function(fun = "norm", n = 1000, robust = TRUE, subdivisions = 100, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Checks consistency of distributions

    # Arguments:
    #   fun - a character string denoting the name of the distribution
    #   n - an integer specifying the number of random variates to be
    #       generated
    #   robust -  a logical flag, should robust estimates be used? By
    #       default \code{TRUE}
    #   subdivisions - an integer specifying the numbers of subdivisions
    #       in integration __ NB: only used in one place, *not* in the other.. hmm

    #   ... - the distributional parameters and passed to integrate()
### FIXME (MM): add  args.integrate = list()   to be passed to integrate(),
### -----       and pass the others to distrib.functions

    # Examples:
    #   .distCheck("norm", mean = 1, sd = 1)
    #   .distCheck("t", df = 4)
    #   .distCheck("exp", rate = 2)
    #   .distCheck("weibull", shape = 1)

    # FUNCTION:

    # Distribution Functions:
    cat("\nDistribution Check for:", fun, "\n ")
    CALL = match.call()
    cat("Call: ")
    cat(paste(deparse(CALL), sep = "\n", collapse = "\n"), "\n", sep = "")
    dfun = match.fun(paste("d", fun, sep = ""))
    pfun = match.fun(paste("p", fun, sep = ""))
    qfun = match.fun(paste("q", fun, sep = ""))
    rfun = match.fun(paste("r", fun, sep = ""))

    # Range:
    xmin = qfun(p = 0.01, ...)
    xmax = qfun(p = 0.99, ...)

    # Check 1 - Normalization:
    NORM = integrate(dfun, lower = -Inf, upper = Inf,
        subdivisions = subdivisions, stop.on.error = FALSE, ...)
    cat("\n1. Normalization Check:\n NORM ")
    print(NORM)
    normCheck = (abs(NORM[[1]]-1) < 0.01)

    # Check 2:
    cat("\n2. [p-pfun(qfun(p))]^2 Check:\n ")
    p = c(0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999)
    P = pfun(qfun(p, ...), ...)
    pP = round(rbind(p, P), 3)
    print(pP)
    RMSE = sd(p-P)
    print(c(RMSE = RMSE))
    rmseCheck = (abs(RMSE) < 0.0001)

    # Check 3:
    cat("\n3. r(", n, ") Check:\n", sep = "")
    r = rfun(n = n, ...)
    if (!robust) {
        SAMPLE.MEAN = mean(r)
        SAMPLE.VAR = var(r)
    } else {
        robustSample = MASS::cov.mcd(r, quantile.used = floor(0.95*n))
        SAMPLE.MEAN = robustSample$center
        SAMPLE.VAR = robustSample$cov[1,1]
    }
    SAMPLE = data.frame(t(c(MEAN = SAMPLE.MEAN, "VAR" = SAMPLE.VAR)),
        row.names = "SAMPLE")
    print(signif(SAMPLE, 3))
    fun1 = function(x, ...) { x * dfun(x, ...) }
    fun2 = function(x, M, ...) { x^2 * dfun(x, ...) }
    MEAN = integrate(fun1, lower = -Inf, upper = Inf,
        subdivisions = 5000, stop.on.error = FALSE,...)
    cat("   X   ")
    print(MEAN)
    VAR = integrate(fun2, lower = -Inf, upper = Inf,
        subdivisions = 5000, stop.on.error = FALSE, ...)
    cat("   X^2 ")
    print(VAR)
    EXACT = data.frame(t(c(MEAN = MEAN[[1]], "VAR" = VAR[[1]] - MEAN[[1]]^2)),
        row.names = "EXACT ")
    print(signif(EXACT, 3))
    meanvarCheck = (abs(SAMPLE.VAR-EXACT$VAR)/EXACT$VAR < 0.1)
    cat("\n")

    # Done:
    ans = list(
        normCheck = normCheck,
        rmseCheck = rmseCheck,
        meanvarCheck = meanvarCheck)

    # Return Value:
    unlist(ans)
}


# ------------------------------------------------------------------------------


.distCheck <- distCheck

    # Keep for older Rmetrics Versions



################################################################################

