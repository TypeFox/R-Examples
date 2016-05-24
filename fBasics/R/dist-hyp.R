
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
#  dhyp                  Returns density for hyperbolic DF
#  phyp                  Returns probability for hyperbolic DF
#  qhyp                  Returns quantiles for hyperbolic DF
#  rhyp                  Returns random variates for hyperbolic DF
################################################################################


dhyp <-
function(x, alpha = 1, beta = 0, delta = 1, mu = 0, pm = c("1", "2", "3", "4"), 
    log = FALSE)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns Hyperbolic Density Function PDF

    # Arguments:
    #   alpha, beta - Shape Parameter, |beta| <= alpha
    #   delta  - Scale Parameter, 0 <= delta
    #   mu - Location Parameter

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
    
    # Settings:
    pm = match.arg(pm)

    # Density:
    ans = NA
    if (pm == 1) ans = .dhyp1(x, alpha, beta, delta, mu)
    if (pm == 2) ans = .dhyp2(x, alpha, beta, delta, mu)
    if (pm == 3) ans = .dhyp3(x, alpha, beta, delta, mu)
    if (pm == 4) ans = .dhyp4(x, alpha, beta, delta, mu)
    
    # Log:
    if (log) ans = log(ans)
    
    # Return value:
    ans
}


# ------------------------------------------------------------------------------


phyp <-
function(q, alpha = 1, beta = 0, delta = 1, mu = 0, pm = c("1", "2", "3", "4"), 
    ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Return cumulative probability of Hyperbolic PDF

    # Arguments:
    #   alpha, beta - Shape Parameter, |beta| <= alpha
    #   delta  - Scale Parameter, 0 <= delta
    #   mu - Location Parameter

    # FUNCTION:

    # Checks:
    if (alpha <= 0) stop("alpha must be greater than zero")
    if (delta <= 0) stop("delta must be greater than zero")
    if (abs(beta) >= alpha) stop("abs value of beta must be less than alpha")
    
    # Settings:
    pm = match.arg(pm)

    # Return Value:
    ans = NA
    if (pm == 1) return(.phyp1(q, alpha, beta, delta, mu, ...))
    if (pm == 2) return(.phyp2(q, alpha, beta, delta, mu, ...))
    if (pm == 3) return(.phyp3(q, alpha, beta, delta, mu, ...))
    if (pm == 4) return(.phyp4(q, alpha, beta, delta, mu, ...))
}


# ------------------------------------------------------------------------------


qhyp <-
function(p, alpha = 1, beta = 0, delta = 1, mu = 0, pm = c("1", "2", "3", "4"), 
    ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns quantiles of Hyperbolic PDF

    # Arguments:
    #   alpha, beta - Shape Parameter, |beta| <= alpha
    #   delta  - Scale Parameter, 0 <= delta
    #   mu - Location Parameter

    # Note:
    #   This procedure will not run under Splus.

    # FUNCTION:

    # Checks:
    if (alpha <= 0) stop("alpha must be greater than zero")
    if (delta <= 0) stop("delta must be greater than zero")
    if (abs(beta) >= alpha) stop("abs value of beta must be less than alpha")
    
    # Settings:
    pm = match.arg(pm)

    # Return Value:
    ans = NA
    if (pm == 1) return(.qhyp1(p, alpha, beta, delta, mu, ...))
    if (pm == 2) return(.qhyp2(p, alpha, beta, delta, mu, ...))
    if (pm == 3) return(.qhyp3(p, alpha, beta, delta, mu, ...))
    if (pm == 4) return(.qhyp4(p, alpha, beta, delta, mu, ...))
}


# ------------------------------------------------------------------------------


rhyp <-
function(n, alpha = 1, beta = 0, delta = 1, mu = 0, pm = c("1", "2", "3", "4"))
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns random deviates of Hyperbolic PDF

    # Arguments:
    #   n - number of random deviates to be generated
    #   alpha, beta - Shape Parameter, |beta| <= alpha
    #   delta  - Scale Parameter, 0 <= delta
    #   mu - Location Parameter

    # Notes:
    #   I have removed my original Fortran program and replaced it by
    #   the dhyperb() function from the HyperbolicDist Package, written
    #   by David Scott, Ai-Wei Lee, Jennifer Tso, Richard Trendall.
    #   License: GPL

    # FUNCTION:

    # Checks:
    if (alpha <= 0) stop("alpha must be greater than zero")
    if (delta <= 0) stop("delta must be greater than zero")
    if (abs(beta) >= alpha) stop("abs value of beta must be less than alpha")
    
    # Settings:
    pm = match.arg(pm)

    # Result:
    ans = NA
    if (pm == 1) ans = .rhyp1(n, alpha, beta, delta, mu)
    if (pm == 2) ans = .rhyp2(n, alpha, beta, delta, mu)
    if (pm == 3) ans = .rhyp3(n, alpha, beta, delta, mu)
    if (pm == 4) ans = .rhyp4(n, alpha, beta, delta, mu)

    # Attributes:
    attr(ans, "control") = c(dist = "hyp", alpha = alpha, beta = beta,
    delta = delta, mu = mu)

    # Return Value:
    ans
}


################################################################################

