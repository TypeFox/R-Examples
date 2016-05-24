
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
#  dght                  Density of the generalized hyperbolic Student-t 
#  pght                  Probability of the GHT
#  dght                  Quantiles of the GHT
#  rght                  Random variates of the GHT
################################################################################


dght <-
function(x, beta = 0.1, delta = 1, mu = 0, nu = 10, log = FALSE)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns density of the generalized hyperbolic Student-t

    # Arguments:

    # Example:
    #   x = (-5):5; dght(x)

    # FUNCTION:

    # Parameters:
    if (length(beta) == 4) {
       nu = beta[4]
       mu = beta[3]
       delta = beta[2]
       beta = beta[1]
    } 
    
    # GH Parameters:
    alpha = abs(beta) + 1e-6
    lambda = -nu/2
    
    # Density:
    ans = dgh(x, alpha, beta, delta, mu, lambda, log = log)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


pght <-
function(q, beta = 0.1, delta = 1, mu = 0, nu = 10)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns probabilities of the generalized hyperbolic Student-t

    # Arguments:
    
    # Example:
    #   q = (-5):5; pght(q)

    # FUNCTION:

    # Parameters:
    if (length(beta) == 4) {
       nu = beta[4]
       mu = beta[3]
       delta = beta[2]
       beta = beta[1]
    } 
    
    # Cumulative Probability:
    ans = NULL
    for (Q in q) {
        Integral = integrate(dght, -Inf, Q, stop.on.error = FALSE,
            beta = beta, delta = delta, mu = mu, nu = nu)
        ans = c(ans, as.numeric(unlist(Integral)[1]) )
    }

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


qght <-
function(p, beta = 0.1, delta = 1, mu = 0, nu = 10)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns quantiles of the generalized hyperbolic Student-t

    # Arguments:
    
    # Example:
    #   p = (1:9)/10; qght(p); round(pght(qght(p)), digits = 4)

    # FUNCTION:

    # Parameters:
    if (length(beta) == 4) {
       nu = beta[4]
       mu = beta[3]
       delta = beta[2]
       beta = beta[1]
    } 
    
    # Internal Functions:
    .froot <- function(x, beta = beta, delta = delta, mu = mu, nu = nu, p) 
    {
        pght(q = x, beta = beta, delta = delta, mu = mu, nu = nu) - p
    }

    # Loop over all p's:
    result = NULL
    for (pp in p) {
        lower = -1
        upper = +1
        counter = 0
        iteration = NA
        while (is.na(iteration)) {
            iteration = .unirootNA(f = .froot, interval = c(lower, 
                upper), beta = beta, delta = delta, mu = mu, 
                nu = nu, p = pp)
            counter = counter + 1
            lower = lower - 2^counter
            upper = upper + 2^counter
        }
        result = c(result, iteration)
    }

    # Return Value:
    ans = result + mu
    ans
}


# ------------------------------------------------------------------------------


rght <- 
function(n, beta = 0.1, delta = 1, mu = 0, nu = 10)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns random Variates of generalized hyperbolic Student-t

    # Arguments:
    
    # Example:
    #   r = rght(10)

    # FUNCTION:
    
    # Parameters:
    if (length(beta) == 4) {
       nu = beta[4]
       mu = beta[3]
       delta = beta[2]
       beta = beta[1]
    } 
    
    # GH Parameters:
    alpha = abs(beta) + 1e-6
    lambda = -nu/2

    # Random Variates:
    x = rgh(n, alpha, beta = beta, delta = delta, mu = mu, lambda) 

    # Return Value:
    x
}


################################################################################
