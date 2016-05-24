
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
#   1999 - 2007, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:               GPD SIMULATION:
#  gpdSim                  Simulates a GPD distributed process
# FUNCTION:               GPD PARAMETER ESTIMATION:
#  'fGPDFIT'               S4 class representation
#  gpdFit                  Fits Parameters of GPD distribution
#  .gpdpwmFit              Fits GPD with probability weighted moments
#  .gpdmleFit              Fits GPD with max log-likelihood approach
#   .gpdLLH                 Computes GPD log-likelihood function
################################################################################


gpdSim <-
    function(model = list(xi = 0.25, mu = 0, beta = 1), n = 1000, seed = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Generates random variates from a GPD distribution

    # FUNCTION:

    # Seed:
    if (is.null(seed)) seed = NA else set.seed(seed)

    # Simulate:
    ans = rgpd(n = n, xi = model$xi, mu = model$mu, beta = model$beta)
    ans = as.ts(ans)

    # Control:
    attr(ans, "control") =
        data.frame(t(unlist(model)), seed = seed, row.names = "")

    # Return Value:
    ans
}


################################################################################


setClass("fGPDFIT",
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


gpdFit <-
    function(x, u = quantile(x, 0.95), type = c("mle", "pwm"),
    information = c("observed", "expected"),
    title = NULL, description = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits a generalized Pareto model to excesses

    # Arguments:

    # Details:
    #   Returns an object of class "fGPDFIT" representing the fit of a
    #   generalized Pareto model to excesses over a high threshold

    # Notes:
    #   This is a wrapper to EVIR's 'gpd' function.

    # FUNCTION:

    # Settings:
    call = match.call()
    type = match.arg(type)
    information = match.arg(information)

    # Check Type and Convert:
    X = x
    xClass = class(x)
    if (xClass == "timeSeries") stopifnot(isUnivariate(x))
    x = as.vector(x)
    N = length(x)

    # Compute Exceedances:
    exceedances = x[x > u]
    Names = as.character((1:N)[x > u])
    exceedances = as.vector(exceedances)
    names(exceedances) = Names

    # Estimate Parameters:
    if (type == "mle") {
        fit = .gpdmleFit(x, u, information)
        fit$llh = fit$fit$value
        fit$convergence = fit$fit$convergence
    } else if (type == "pwm") {
        fit = .gpdpwmFit(x, u)
        fit$llh = NA
        fit$convergence = NA
    }
    fit$p.less.thresh = fit$prob = 1 - length(x[x > u]) / length(x)
    fit$threshold = u
    fit$data = x

    # Compute Residuals:
    xi = fit$par.ests["xi"]
    beta = fit$par.ests["beta"]
    residuals = log(1 + (xi * (as.vector(exceedances)-u))/beta)/xi

    # Add title and description:
    if (is.null(title)) title = "GPD Parameter Estimation"
    if (is.null(description)) description = description()

    # Compose Parameter List:
    parameter = list(u = u, type = type)
    if (information == "mle") parameter$information = information

    # Return Value:
    new("fGPDFIT",
        call = call,
        method = c("gpd", type),
        parameter = parameter,
        data = list(x = X, exceedances = exceedances),
        fit = fit,
        residuals = residuals,
        title = title,
        description = description)
}


# ------------------------------------------------------------------------------


.gpdmleFit <-
    function (x, u, information = c("observed", "expected"), ...)
{
    # A Copy from Evir

    data = x
    threshold = u

    exceedances <- data[data > threshold]
    excess <- exceedances - threshold
    Nu <- length(excess)
    xbar <- mean(excess)

    s2 <- var(excess)
    xi0 <- -0.5 * (((xbar * xbar)/s2) - 1)
    beta0 <- 0.5 * xbar * (((xbar * xbar)/s2) + 1)
    theta <- c(xi0, beta0)

    negloglik <- function(theta, tmp)
    {
        xi <- theta[1]
        beta <- theta[2]
        cond1 <- beta <= 0
        cond2 <- (xi <= 0) && (max(tmp) > (-beta/xi))
        if (cond1 || cond2) {
            f <- 1e+06
        } else {
            y <- logb(1 + (xi * tmp)/beta)
            y <- y/xi
            f <- length(tmp) * logb(beta) + (1 + xi) * sum(y)
        }
        f
    }

    fit <- optim(theta, negloglik, hessian = TRUE, ..., tmp = excess)
    names(fit$par) = c("xi", "beta")

    if (fit$convergence) warning("Optimization may not have been succeeded.")
    par.ests <- fit$par
    converged <- fit$convergence
    nllh.final <- fit$value

    information <- match.arg(information)
    if (information == "observed")
        varcov <- solve(fit$hessian)
    if (information == "expected") {
        one <- (1 + par.ests[1])^2/Nu
        two <- (2 * (1 + par.ests[1]) * par.ests[2]^2)/Nu
        cov <- -((1 + par.ests[1]) * par.ests[2])/Nu
        varcov <- matrix(c(one, cov, cov, two), 2)
    }

    par.ses <- sqrt(diag(varcov))
    names(par.ses) = c("xi", "beta")

    list(par.ests = par.ests, par.ses = par.ses, fit = fit, varcov = varcov)
}


# ------------------------------------------------------------------------------


.gpdmleFitCheck =
    function(x, u = quantile(x, 0.95),
    information = c("observed", "expected"), ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits GPD with max log-likelihood approach

    # FUNCTION:

    x = as.vector(x)
    excess = x[x > u] - u
    theta = .gpdpwmFit(x = x, u = u)$par.ests

    # Parameter Estimation:
    fit = optim(theta, .gpdLLH, hessian = TRUE, excess = excess)
    names(fit$par.ests) = c("xi", "beta")

    # Error Estimates:
    if (information[1] == "observed") {
        varcov = solve(fit$hessian)
    }
    if (information[1] == "expected") {
        Nu = length(excess)
        one = (1 + fit$par[1])^2/Nu
        two = (2 * (1 + fit$par[1]) * fit$par[2]^2)/Nu
        cov = -((1 + fit$par[1]) * fit$par[2])/Nu
        varcov = matrix(c(one, cov, cov, two), 2)
    }
    par.ses = sqrt(diag(varcov))
    names(par.ses) = c("xi", "beta")

    # Return Value:
    list(par.ests = fit$par, par.ses = par.ses, fit = fit, varcov = varcov)
}


# ------------------------------------------------------------------------------


.gpdpwmFit <-
    function (x, u)
{
    # A Copy from Evir

    data = x
    threshold = u

    data <- as.numeric(data)
    n <- length(data)
    exceedances <- data[data > threshold]
    excess <- exceedances - threshold
    Nu <- length(excess)
    xbar <- mean(excess)

    a0 <- xbar
    gamma <- -0.35
    delta <- 0
    pvec <- ((1:Nu) + gamma)/(Nu + delta)
    a1 <- mean(sort(excess) * (1 - pvec))
    xi <- 2 - a0/(a0 - 2 * a1)
    beta <- (2 * a0 * a1)/(a0 - 2 * a1)
    par.ests <- c(xi, beta)
    names(par.ests) = c("xi", "beta")

    denom <- Nu * (1 - 2 * xi) * (3 - 2 * xi)
    if (xi > 0.5) {
        denom <- NA
        warning("Asymptotic Standard Errors not available for PWM when xi>0.5.")
    }
    one <- (1 - xi) * (1 - xi + 2 * xi^2) * (2 - xi)^2
    two <- (7 - 18 * xi + 11 * xi^2 - 2 * xi^3) * beta^2
    cov <- beta * (2 - xi) * (2 - 6 * xi + 7 * xi^2 - 2 * xi^3)
    varcov <- matrix(c(one, cov, cov, two), 2)/denom
    information <- "expected"
    converged <- NA
    nllh.final <- NA
    par.ses <- sqrt(diag(varcov))
    names(par.ses) = c("xi", "beta")

    list(par.ests = par.ests, par.ses = par.ses, fit = NA, varcov = NA)
}



# ------------------------------------------------------------------------------


.gpdpwmFitCheck <-
    function(x, u = quantile(x, 0.95))
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits GPD with probability weighted moments

    # FUNCTION:

    # PWM Fit:
    x = as.vector(x)
    excess = x[x > u] - u
    Nu = length(excess)
    a0 = mean(excess)
    pvec = ((1:Nu) - 0.35)/Nu
    a1 = mean(sort(excess) * (1 - pvec))
    xi = 2 - a0/(a0 - 2 * a1)
    beta = (2 * a0 * a1)/(a0 - 2 * a1)
    par.ests = c(xi = xi, beta = beta)
    names(par.ests) = c("xi", "beta")
    par.ses = c(xi = NA, beta = NA)
    names(par.ses) = c("xi", "beta")

    # Return Value:
    list(par.ests = par.ests, par.ses = par.ses, fit = NA, varcov = NA)
}


################################################################################


.gpdLLH =
function(theta, excess)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes GPD log-likelihood function

    # FUNCTION:

    # LLH:
    xi = theta[1]
    beta = theta[2]
    cond = (beta <= 0) || ((xi <= 0) && (max(excess) > (-beta/xi)))
    if (cond) {
        func = NA
    } else {
        y = log(1+(xi*excess)/beta) / xi
        func = length(excess) * log(beta) + (1+xi)*sum(y)
    }

    # Return Value:
    func
}


################################################################################

