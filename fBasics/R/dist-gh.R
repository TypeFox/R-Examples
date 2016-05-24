
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
#  dgh                   Returns density for generalized hyperbolic DF
#  pgh                   Returns probability for generalized hyperbolic DF
#  qgh                   Returns quantiles for generalized hyperbolic DF
#  rgh                   Returns random variates for generalized hyperbolic DF
################################################################################


dgh <-
function(x, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = -1/2, log = FALSE)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns density for the generalized hyperbolic distribution

    # FUNCTION:

    # Parameters:
    if (length(alpha) == 4) {
       mu = alpha[4]
       delta = alpha[3]
       beta = alpha[2]
       alpha = alpha[1]
    } 
    
    # Checks:
    if (alpha <= 0) stop("alpha must be greater than zero")
    if (delta <= 0) stop("delta must be greater than zero")
    if (abs(beta) >= alpha) stop("abs value of beta must be less than alpha")

    # Density:
    arg = delta*sqrt(alpha^2-beta^2)
    a = (lambda/2)*log(alpha^2-beta^2) - (
        log(sqrt(2*pi)) + (lambda-0.5)*log(alpha) + lambda*log(delta) +
        log(besselK(arg, lambda, expon.scaled = TRUE)) - arg )
    f = ((lambda-0.5)/2)*log(delta^2+(x - mu)^2)

    # Use exponential scaled form to prevent from overflows:
    arg = alpha * sqrt(delta^2+(x-mu)^2)
    k = log(besselK(arg, lambda-0.5, expon.scaled = TRUE)) - arg
    e = beta*(x-mu)

    # Put all together:
    ans = a + f + k + e
    if(!log) ans = exp(ans)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


pgh <-
function(q, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = -1/2)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns probability for the generalized hyperbolic distribution

    # FUNCTION:

    # Checks:
    if (alpha <= 0) stop("alpha must be greater than zero")
    if (delta <= 0) stop("delta must be greater than zero")
    if (abs(beta) >= alpha) stop("abs value of beta must be less than alpha")

    # Probability:
    ans = NULL
    for (Q in q) {
        Integral = integrate(dgh, -Inf, Q, stop.on.error = FALSE,
            alpha = alpha, beta = beta, delta = delta, mu = mu,
            lambda = lambda)
        ans = c(ans, as.numeric(unlist(Integral)[1]))
    }

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


qgh <-
function(p, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = -1/2)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns quantiles for the generalized hyperbolic distribution

    # FUNCTION:

    # Checks:
    if (alpha <= 0) stop("alpha must be greater than zero")
    if (delta <= 0) stop("delta must be greater than zero")
    if (abs(beta) >= alpha) stop("abs value of beta must be less than alpha")

    # Internal Function:
    .froot <- function(x, alpha, beta, delta, mu, lambda, p)
    {
        pgh(q = x, alpha = alpha, beta = beta, delta = delta,
            mu = mu, lambda = lambda) - p
    }

    # Quantiles:
    result = NULL
    for (pp in p) {
        lower = -1
        upper = +1
        counter = 0
        iteration = NA
        while (is.na(iteration)) {
            iteration = .unirootNA(f = .froot, interval = c(lower,
                upper), alpha = alpha, beta = beta, delta = delta,
                mu = mu, lambda = lambda, p = pp)
            counter = counter + 1
            lower = lower - 2^counter
            upper = upper + 2^counter
        }
        result = c(result, iteration)
    }

    # Return Value:
    result
}


# ------------------------------------------------------------------------------


rgh <-
function(n, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = -1/2)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns random variates for the generalized hyperbolic distribution

    # FUNCTION:

    # Checks:
    if (alpha <= 0) stop("alpha must be greater than zero")
    if (delta <= 0) stop("delta must be greater than zero")
    if (abs(beta) >= alpha) stop("abs value of beta must be less than alpha")

    # Settings:
    theta = c(lambda, alpha, beta, delta, mu)

    # Random Numbers:
    ans = .rghyp(n, theta)

    # Attributes:
    attr(ans, "control") = c(dist = "gh", alpha = alpha, beta = beta,
    delta = delta, mu = mu, lambda = lambda)

    # Return Value:
    ans
}


################################################################################

