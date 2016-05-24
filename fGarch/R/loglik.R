
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
# FUNCTION:               DESCRIPTION:
#  .garchLLH               Computes log-likelihood function
#  .garchOptimizeLLH       Opimizes log-likelihood function
################################################################################


.garchLLH <-
    function(params, trace = TRUE, fGarchEnv = FALSE)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Compute Log-Likelihood Function

    # Arguments:
    #   params - a named numeric vector with the model parameters
    #       to be optimized

    # Value:
    #   Returns the value of the max log-likelihood function.

    # Note:
    #   The variables '.series' and '.params' must be global available

    # FUNCTION:

    # DEBUG:
    DEBUG = FALSE
    if (DEBUG) print("Entering Function .garchLLH")

    # Get Global Variables:
    .series <- .getfGarchEnv(".series")
    .params <- .getfGarchEnv(".params")
    .garchDist <- .getfGarchEnv(".garchDist")
    .llh <- .getfGarchEnv(".llh")

    # How to calculate the LLH Function?
    if (DEBUG) print(.params$control$llh)

    if(.params$control$llh == "internal")  {

        if (DEBUG) print("internal")
        return(.aparchLLH.internal(params, trace = trace, fGarchEnv = fGarchEnv))

    } else if (.params$control$llh == "filter")  {

        if (DEBUG) print("filter")
        return(.aparchLLH.filter(params, trace = trace, fGarchEnv = fGarchEnv))

    } else if (.params$control$llh == "testing")  {

        if (DEBUG) print("testing")
        return(.aparchLLH.testing(params, trace = trace, fGarchEnv = fGarchEnv))

    } else {

        stop("LLH is neither internal, testing, nor filter!")
    }

}


# ------------------------------------------------------------------------------


.garchOptimizeLLH <-
function(hessian = hessian, robust.cvar, trace)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Opimizes the Log-Likelihood Function

    # Arguments:
    #   hessian - the Hessian matrix
    #   robust.cvar - a logical
    #   trace - a logical

    # FUNCTION:

    # DEBUG:
    DEBUG = FALSE
    if (DEBUG) print("Entering Function .garchOptimizeLLH")

    # get global variables
    .series <- .getfGarchEnv(".series")
    .params <- .getfGarchEnv(".params")

    # Initialization:
    INDEX = .params$index

    # Algorithm:
    algorithm = .params$control$algorithm[1]
    TOL1 = .params$control$tol1
    TOL2 = .params$control$tol2
    if(trace) {
        cat("\n\n--- START OF TRACE ---")
        cat("\nSelected Algorithm:", algorithm, "\n")
    }

    # First Method:
    # Two Step Apparoach > Trust Region + Nelder-Mead Simplex
    if(algorithm == "nlminb" | algorithm == "nlminb+nm") {
        fit <- .garchRnlminb(.params, .series, .garchLLH, trace)
        .params$llh = fit$llh
        .params$params[INDEX] = fit$par
        .setfGarchEnv(.params = .params)
    }
    if(algorithm == "nlminb+nm") {
        fit <- .garchRnm(.params, .series, .garchLLH, trace)
            .params$llh = fit$llh
            .params$params[INDEX] = fit$par
            .setfGarchEnv(.params = .params)
    }

    # Second Method:
    # Two Step Approach > BFGS + Nelder-Mead Simplex
    if(algorithm == "lbfgsb" | algorithm == "lbfgsb+nm") {
        fit <- .garchRlbfgsb(.params, .series, .garchLLH, trace)
        .params$llh = fit$llh
        .params$params[INDEX] = fit$par
        .setfGarchEnv(.params = .params)
    }
    if(algorithm == "lbfgsb+nm") {
        fit <- .garchRnm(.params, .series, .garchLLH, trace)
        .params$llh = fit$llh
        .params$params[INDEX] = fit$par
        .setfGarchEnv(.params = .params)
    }

    # Save parameters:
    .params$llh = fit$llh
    .params$params[INDEX] = fit$par
    .setfGarchEnv(.params = .params)

    # Compute the Hessian:
    if (hessian == "ropt") {
        fit$hessian <- - .garchRoptimhess(par = fit$par, .params = .params,
            .series = .series)
        titleHessian = "R-optimhess"
    } else if (hessian == "rcd") {
        fit$hessian <- - .garchRCDAHessian(par = fit$par, .params = .params,
            .series = .series)
        titleHessian = "Central"
    } else if (hessian == "rts") {
        fit$hessian <- - .garchTSHessian(par = fit$par, .params = .params,
            .series = .series)
        titleHessian = "Two Sided"
    }

    # Rescale Parameters:
    if (.params$control$xscale) {
        .series$x <- .series$x * .series$scale
        if (.params$include["mu"])
            fit$coef["mu"] <- fit$par["mu"] <- .params$params["mu"] <-
                .params$params["mu"]*.series$scale
        if (.params$include["omega"])
            fit$coef["omega"] <- fit$par["omega"] <- .params$params["omega"] <-
                .params$params["omega"]*.series$scale^(.params$params["delta"])
        # save changes
        .setfGarchEnv(.params = .params)
        .setfGarchEnv(.series = .series)
    }

    # Rescale Hessian Matrix:
    if (.params$control$xscale) {
        if (.params$include["mu"]) {
            fit$hessian[,"mu"] <- fit$hessian[,"mu"] /  .series$scale
            fit$hessian["mu",] <- fit$hessian["mu",] /  .series$scale
        }
        if (.params$include["omega"]) {
            fit$hessian[,"omega"] <-
                fit$hessian[,"omega"] / .series$scale^(.params$params["delta"])
            fit$hessian["omega",] <-
                fit$hessian["omega",] / .series$scale^(.params$params["delta"])
        }
    }

    # Recalculate llh, h, z with Rescaled Parameters:
    .llh <- fit$llh <- fit$value <-
        .garchLLH(fit$par, trace = FALSE, fGarchEnv = TRUE)
    .series <- .getfGarchEnv(".series")

    # Compute the Gradient:
    # YC: needs to be after the calculation of h, z !
    if (robust.cvar)
        fit$gradient <- - .garchRCDAGradient(
            par = fit$par, .params = .params, .series = .series)

    # Compute Information Criterion Statistics:
    N = length(.series$x)
    NPAR = length(fit$par)
    fit$ics = c(
        AIC  = c((2*fit$value)/N + 2 * NPAR/N),
        BIC  = (2*fit$value)/N + NPAR * log(N)/N,
        SIC  = (2*fit$value)/N + log((N+2*NPAR)/N),
        HQIC = (2*fit$value)/N + (2*NPAR*log(log(N)))/N )
    names(fit$ics) <- c("AIC", "BIC", "SIC", "HQIC")

    # Print LLH if we trace:
    if(trace) {
        cat("\nFinal Estimate of the Negative LLH:\n")
        cat(" LLH: ", .llh, "   norm LLH: ", .llh/N, "\n")
        print(fit$par)
    }

    # Print Hessian Matrix if we trace:
    if(trace) {
        cat("\n", titleHessian, " Difference Approximated Hessian Matrix:\n",
            sep = "")
        print(fit$hessian)
        cat("\n--- END OF TRACE ---\n\n")
    }

    # Return Value:
    if (DEBUG) print("Entering Function .garchOptimizeLLH")
    fit
}


################################################################################

