
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

# Copyrights (C)
# for this R-port:
#   1999 - 2009, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:             GEV PARAMETER ESTIMATION:
#  'fGEVFIT'             S4 class representation
#  gevFit                Fits parameters of GEV distribution
#  gumbelFit             Fits parameters of Gumbel distribution
# FUNCTION:             FOR INTERNAL USE:
#  .gumpwmFit            Fits Gumbel with probability weighted moments
#  .gevpwmFit            Fits GEV with probability weighted moments
#  .gummleFit            Fits Gumbel with max log-likelihood approach
#   .gumLLH               Computes Gumbel log-likelihood function
#  .gevmleFit            Fits GEV with max log-likelihood approach
#   .gevLLH               Computes GEV log-likelihood function
################################################################################


setClass("fGEVFIT",
    representation(
        call = "call",
        method = "character",
        parameter = "list",
        data = "list",
        fit = "list",
        residuals = "numeric",
        title = "character",
        description = "character"
    )
)


# ------------------------------------------------------------------------------


gevFit =
function(x, block = 1, type = c("mle", "pwm"),
title = NULL, description = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits the parameters of GEV distribution

    # Arguments:
    #   x - an object of class timeSeries
    #   block - an integer value, the block size
    #   type - a character string, which type of method should be used,
    #       max log-likelihood estimation, "mle", or partial weighted
    #       moments estimation, "pwm".

    # Examples:
    #   gevFit(gevSim())
    #   par(mfrow = c(2,2)); summary(gevFit(gevSim()))

    # FUNCTION:

    # Match Call:
    call = match.call()

    # Match Arguments:
    type = match.arg(type)

    # Fit:
    ans = .gevFit(x = x, block = block, type = type, gumbel = FALSE,
        title = title, description = description, ...)
    ans@call = call

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


gumbelFit =
function(x, block = 1, type = c("mle", "pwm"),
title = NULL, description = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits the parameters of Gumbel distribution

    # Arguments:
    #   x - an object of class timeSeries
    #   block - an integer value, the block size
    #   type - a character string, which type of method should be used,
    #       max log-likelihood estimation, "mle", or partial weighted
    #       moments estimation, "pwm".

    # Examples:
    #   gumbelFit(gumbelSim())
    #   par(mfrow = c(2,2)); summary(gumbelFit(gumbelSim()))

    # FUNCTION:

    # Match Call:
    call = match.call()

    # Match Arguments:
    type = match.arg(type)

    # Fit:
    ans = .gevFit(x = x, block = block, type = type, gumbel = TRUE,
        title = title, description = description, ...)
    ans@call =  call

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.gevFit =
function(x, block = 1, type = c("mle", "pwm"),
gumbel = FALSE, title = NULL, description = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits parameters to a GEV distribution

    # Arguments:
    #   x - a numeric vector of Block Maxima

    # Examples:
    #   fit = gevFit(gevSim(), type = "mle", gumbel = FALSE); print(fit)
    #   fit = gevFit(gevSim(), type = "pwm", gumbel = FALSE); print(fit)
    #   fit = gevFit(gevSim(), type = "mle", gumbel = TRUE); print(fit)
    #   fit = gevFit(gevSim(), type = "pwm", gumbel = TRUE); print(fit)
    #   x = rnorm(500);  block = 20; type="mle"; gumbel=FALSE
    #   x = as.ts(rnorm(500));  block = 20; type = "mle"; gumbel = FALSE
    #   x = dummyDailySeries(rnorm(500)); block = 20; type = "mle"; gumbel=FALSE

    # Note:
    #   Argument named "method is already used for the selection
    #   of the MLE optimization algorithm, therfore we use here
    #   "type".

    # FUNCTION:

    # Match Call:
    call = match.call()

    # Match Arguments:
    type = match.arg(type)

    # Check Type and Convert:
    X = x
    xClass = class(x)
    x = as.timeSeries(x)
    stopifnot(isUnivariate(x))

    # Block Maxima:
    if (is.numeric(block)) {
        if (block == 1) {
            blockmaxima = x
            Names = paste(1:dim(blockmaxima)[1])
        } else {
            blockmaxima = blockMaxima(x, block, doplot = FALSE)
            Names = blockmaxima@recordIDs[, 3]
        }
    } else {
        blockmaxima = blockMaxima(x, block, doplot = FALSE)
        Names = rownames(series(blockmaxima))
    }

    if (xClass == "numeric") {
        blockmaxima = as.vector(blockmaxima)
        names(blockmaxima) = Names
    }
    if (xClass == "ts") {
        blockmaxima = as.ts(blockmaxima)
        names(blockmaxima) = Names
    }
    x = as.vector(blockmaxima)

    # Estimate Parameters:
    if (gumbel) {
        # GUMBEL: Add Call and Type
        if (length(type) > 1) type = type
        # Probability Weighted Moment Estimation:
        if (type == "pwm") {
            fit = .gumpwmFit(data = x, ...)
        }
        # Maximum Log Likelihood Estimation:
        # Use Alexander McNeils EVIS from evir Package ...
        if (type == "mle") {
            fit = .gummleFit(data = x, ...)
        }
    } else {
        # GEV: Add Call and Type
        if (length(type) > 1) type = type
        # Probability Weighted Moment Estimation:
        if (type == "pwm") {
            fit = .gevpwmFit(data = x, ...)
        }
        # Maximum Log Likelihood Estimation:
        # Use Alexander McNeils EVIS from evir Package ...
        if (type == "mle") {
            fit = .gevmleFit(data = x, ...)
        }
    }
    class(fit) = "list"

    # Compute Residuals:
    if (gumbel) {
        # GUMBEL:
        xi = 0
        beta = fit$par.ests["beta"]
        mu = fit$par.ests["mu"]
        residuals = exp( - exp( - (x - mu)/beta))
    } else {
        # GEV:
        xi = fit$par.ests["xi"]
        beta = fit$par.ests["beta"]
        mu = fit$par.ests["mu"]
        residuals = (1 + (xi * (x - mu))/beta)^(-1/xi)
    }

    # Make Unique:
    fit$llh = fit$nllh.final

    # Add title and description:
    if (is.null(title)) {
        if (gumbel) {
            title = "Gumbel Parameter Estimation"
        } else {
            title = "GEV Parameter Estimation"
        }
    }
    if (is.null(description)) {
        description = as.character(date())
    }

    # Add Counts to x:

    # Return Value:
    new("fGEVFIT",
        call = match.call(),
        method = c(if (gumbel) "gum" else "gev", type[1]),
        parameter = list(block = block, type = type[1], gumbel = gumbel),
        data = list(x = X, blockmaxima = blockmaxima),
        fit = fit,
        residuals = residuals,
        title = title,
        description = description)
}


# ------------------------------------------------------------------------------


.gumpwmFit =
function(data, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:

    # Arguments:

    # FUNCTION:

    # "Probability Weighted Moment" method.
    data = as.numeric(data)
    n = length(data)

    # Sample Moments:
    x = rev(sort(data))
    lambda = c(mean(x), 0)
    for (i in 1:n) {
        weight = (n-i)/(n-1)/n
        lambda[2] = lambda[2] + weight*x[i]
    }

    # Calculate Parameters:
    xi = 0
    beta = lambda[2]/log(2)
    mu = lambda[1] - 0.5772*beta

    # Output:
    fit = list(
        n = n,
        data = data,
        par.ests = c(mu = mu, beta = beta),
        par.ses = c(mu = NA, beta = NA),
        varcov = matrix(rep(NA, 4), 2, 2),
        converged = NA,
        nllh.final = NA,
        call = match.call(),
        selected = "pwm")
    class(fit) = "gev" # not gumbel!

    # Return Value:
    fit
}


# ------------------------------------------------------------------------------


.gevpwmFit =
function(data, block = NA, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:

    # Arguments:

    # FUNCTION:

    # Probability Weighted Moment method.
    data = as.numeric(data)
    n = length(data)

    # Internal Function:
    y = function(x, w0, w1, w2) {
        (3^x-1)/(2^x-1) - (3*w2 - w0)/(2*w1 - w0)
    }

    # Moments:
    nmom = 3
    x = rev(sort(data))
    moments = rep(0, nmom)
    moments[1] = mean(x)
    n = length(x)
    for (i in 1:n) {
        weight = 1/n
        for (j in 2:nmom) {
            weight = weight*(n-i-j+2)/(n-j+1)
            moments[j] = moments[j] + weight*x[i]
        }
    }
    w0 = moments[1]
    w1 = moments[2]
    w2 = moments[3]

    # Parameters:
    xi = uniroot(f = y, interval = c(-5,+5),
        w0 = w0, w1 = w1, w2 = w2)$root
    beta = (2*w1-w0)*xi / gamma(1-xi) / (2^xi-1)
    mu = w0 + beta*(1-gamma(1-xi))/xi

    # Output:
    fit = list(
        n = n,
        data = data,
        par.ests = c(xi = xi, mu = mu, beta = beta),
        par.ses = c(xi = NA, mu = NA, beta = NA),
        varcov = matrix(rep(NA, 9), 3, 3),
        converged = NA,
        nllh.final = NA,
        call = match.call(),
        selected = "pwm")
    class(fit) = "gev"

    # Return Value:
    fit
}


# ------------------------------------------------------------------------------


.gummleFit =
function(data, block = NA, ...)
{
    # A copy from evir

    # Description:

    # Arguments:

    # FUNCTION:

    # Data:
    data = as.numeric(data)
    n = length(data)

    # Generate EVIR Start Values:
    # beta0 = sqrt(6 * var(data))/pi
    # mu0 = mean(data) - 0.57722 * beta0
    # theta = c(mu = mu0, beta = beta0)
    # We use PWM Start Values:
    theta = .gumpwmFit(data)$par.ests

    # Fit:
    fit = optim(theta, .gumLLH, hessian = TRUE, ..., tmp = data)
    if( fit$convergence) warning("optimization may not have succeeded")
    par.ests = fit$par
    varcov = solve(fit$hessian)
    par.ses = sqrt(diag(varcov))

    # Result:
    ans = list(
        n = n,
        data = data,
        par.ests = par.ests,
        par.ses = par.ses,
        varcov = varcov,
        converged = fit$convergence,
        nllh.final = fit$value)
    class(ans) = "gev"

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.gumLLH =
function(theta, tmp)
{
    # A copy from evir

    # Description:

    # Arguments:

    # FUNCTION:

    # Gumbel Log-Likelihood:
    y = (tmp - theta[1])/theta[2]
    if(theta[2] < 0) {
        ans = 1.0e+6
    } else {
        term1 = length(tmp) * logb(theta[2])
        term2 = sum(y)
        term3 = sum(exp( - y))
        ans = term1 + term2 + term3
    }

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.gevmleFit =
function(data, block = NA, ...)
{
    # A copy from evir

    # Description:

    # Arguments:

    # FUNCTION:

    # Data:
    data = as.numeric(data)
    n = length(data)

    # EVIR Start Values:
    beta0 = sqrt(6 * var(data))/pi
    mu0 = mean(data) - 0.57722 * beta0
    xi0 = 0.1

    # We use PWM Start Values:
    theta = .gevpwmFit(data)$par.ests

    # Fit:
    fit = optim(theta, .gevLLH, hessian = TRUE, ..., tmp = data)
    if (fit$convergence) warning("optimization may not have succeeded")
    par.ests = fit$par
    varcov = solve(fit$hessian)
    par.ses = sqrt(diag(varcov))

    # Result:
    ans = list(
        n = n,
        data = data,
        par.ests = par.ests,
        par.ses = par.ses,
        varcov = varcov,
        converged = fit$convergence,
        nllh.final = fit$value)
    class(ans) = "gev"

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.gevLLH =
function(theta, tmp)
{
    # A copy from evir

    # Description:
    #   Computes log-likelihood for GEV distribution

    # Arguments:

    # FUNCTION:

    # GEV Log-likelihood:
    y = 1 + (theta[1] * (tmp - theta[2]))/theta[3]
    if((theta[3] < 0) || (min(y) < 0)) {
        ans = 1e+06
    } else {
        term1 = length(tmp) * logb(theta[3])
        term2 = sum((1 + 1/theta[1]) * logb(y))
        term3 = sum(y^(-1/theta[1]))
        ans = term1 + term2 + term3
    }

    # Return Value:
    ans
}


################################################################################

