
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

# Copyrights (C)
# for this R-port:
#   1999 - 2008, Diethelm Wuertz, Rmetrics Foundation, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:               SOLVER:
#  .garchRnlminb           R coded solver nlmin
#  .garchRlbgfsb           R coded solver optim using method lbgfsb
#  .garchRnm               R coded solver nm as hybrid addon
################################################################################


.garchRnlminb <-
    function(.params, .series, .garchLLH, trace)
{

    # A function implemented by Diethelm Wuertz

    # Description:
    #   Port3 nlminb R-code Solver

    # Arguments:

    # FUNCTION:

    # Port3 nlminb R-code Solver:
    if(trace) cat("\nR coded nlminb Solver: \n\n")

    # Scale Function and Parameters:
    INDEX = .params$index
    parscale = rep(1, length = length(INDEX))
    names(parscale) = names(.params$params[INDEX])
    parscale["omega"] = var(.series$x)^(.params$delta/2)
    parscale["mu"] = abs(mean(.series$x))

    # Control:
    TOL1 = .params$control$tol1

    # Fit Parameters - par | objective:
    fit <- nlminb(
        start = .params$params[INDEX],
        objective = .garchLLH,
        lower = .params$U[INDEX],
        upper = .params$V[INDEX],
        scale = 1/parscale,
        control = list(
            eval.max = 2000,
            iter.max = 1500,
            rel.tol = 1.0e-14 * TOL1,
            x.tol = 1.0e-14 * TOL1,
            trace = as.integer(trace)),
        fGarchEnv = FALSE) # to speed up .garchLLH
    fit$value <- fit.llh <- fit$objective
    names(fit$par) = names(.params$params[INDEX])

    # make sure to save new h and z (important to speed up code)
    .garchLLH(fit$par, trace = FALSE, fGarchEnv = TRUE)

    # Result:
    fit$coef <- fit$par
    fit$llh <- fit$objective

    # Return Value:
    fit
}


# ------------------------------------------------------------------------------


.garchRlbfgsb <-
    function(.params, .series, .garchLLH, trace)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   optim[L-BFGS-B] Solver

    # Arguments:

    # FUNCTION:

    # optim[L-BFGS-B] Solver:
    if(trace) cat("\nR coded optim[L-BFGS-B] Solver: \n\n")

    # Scale Function and Parameters:
    INDEX = .params$index
    parscale = rep(1, length = length(INDEX))
    names(parscale) = names(.params$params[INDEX])
    parscale["omega"] = var(.series$x)^((.params$params["delta"])/2)

    # Control:
    TOL1 = .params$control$tol1

    # Fit Parameters - par, value:
    fit <- optim(
        par = .params$params[INDEX],
        fn = .garchLLH,
        lower = .params$U[INDEX],
        upper = .params$V[INDEX],
        method = "L-BFGS-B",
        control = list(
            parscale = parscale,
            lmm = 20,
            pgtol = 1.0e-11 * TOL1,
            factr = 1.0 * TOL1,
            trace = as.integer(trace)),
        fGarchEnv = FALSE)  # to speed up .garchLLH
    names(fit$par) = names(.params$params[INDEX])

    # make sure to save new h and z (important to speed up code)
    .garchLLH(fit$par, trace = FALSE, fGarchEnv = TRUE)

    # Print Hessian Matrix:
    # print(fit$hessian)

    # Add to Result:
    fit$coef = fit$par
    fit$llh = fit$value

    # Return Value:
    fit
}


# ------------------------------------------------------------------------------


.garchRnm <-
    function(.params, .series, .garchLLH, trace)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Nelder-Mead as Hybrid Solver

    # Arguments:

    # FUNCTION:

    # Nelder-Mead as Hybrid Solver:
    if(trace) cat("\nR coded Nelder-Mead Hybrid Solver: \n\n")

    # Scale Function and Parameters:
    INDEX = .params$index
    fnscale = abs(.params$llh)
    parscale = abs(.params$params[INDEX])

    # Control:
    TOL2 = .params$control$tol2

    # Fit Parameters:
    fit = optim(
        par = .params$params[INDEX],
        fn = .garchLLH,
        method = "Nelder-Mead",
        control = list(
            ndeps = rep(1e-14*TOL2, length = length(INDEX)),
            maxit = 10000,
            reltol = 1.0e-11 * TOL2,
            fnscale = fnscale,
            parscale = c(1, abs((.params$params[INDEX])[-1])),
            trace = as.integer(trace)),
            hessian = TRUE,
            fGarchEnv = FALSE)  # to speed up .garchLLH
    names(fit$par) = names(.params$params[INDEX])

    # Make sure to save new h and z: (important to speed up code)
    .garchLLH(fit$par, trace = FALSE, fGarchEnv = TRUE)

    # Result:
    fit$coef = fit$par
    fit$llh = fit$value

    # Return Value:
    fit
}

################################################################################

