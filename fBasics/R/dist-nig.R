
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
#  dnig                  Returns density for inverse Gaussian DF
#  pnig                  Returns probability for for inverse Gaussian DF
#  qnig                  Returns quantiles for for inverse Gaussian DF
#  rnig                  Returns random variates for inverse Gaussian DF
# FUNCTION:             DESCRIPTION:
#  .pnigC                Fast C implementation
#  .qnigC                Fast C implementation
################################################################################

.dnig <- function(x, alpha, beta, delta, mu, log = FALSE)
{
    x. <- x - mu
    log.a <- delta*sqrt(alpha^2-beta^2) + log(delta*alpha/pi)
    Sqrt <- sqrt(delta^2 + x.^2)
    log.K1 <- log(besselK(alpha * Sqrt, 1, expon.scaled = TRUE)) - alpha*Sqrt
    dnig <- log.a -log(Sqrt) + log.K1 + beta*x.

    if(log) dnig else exp(dnig)
}

dnig <-
function(x, alpha = 1, beta = 0, delta = 1, mu = 0, log = FALSE)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns density for inverse Gaussian DF

    # Example:
    #   x = rnorm(1000); u = dgh(x, 1.1, 0.2, 0.8, 0.4, -0.5)
    #   v = dnig(rnorm(x), 1.1, 0.2, 0.8, 0.4); u -v

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

    .dnig(x, alpha=alpha, beta=beta, delta=delta, mu=mu, log=log)
}


# ------------------------------------------------------------------------------


pnig <-
function(q, alpha = 1, beta = 0, delta = 1, mu = 0)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns probability for for inverse Gaussian DF

    # Function:

    # Probability:
    #   pgh(q = q, alpha = alpha, beta = beta, delta = delta, mu = mu,
    #   lambda = -0.5)

    if (alpha <= 0)
        stop("alpha must be greater than zero")
    if (delta <= 0)
        stop("delta must be greater than zero")
    if (abs(beta) >= alpha)
        stop("abs value of beta must be less than alpha")
    ans <- q
    for(i in seq_along(q)) {
        ans[i] <- integrate(.dnig, -Inf, q[i],
                            stop.on.error = FALSE,
                            alpha = alpha, beta = beta,
                            delta = delta, mu = mu)$value
    }
    ans
}


# ------------------------------------------------------------------------------


qnig <-
function(p, alpha = 1, beta = 0, delta = 1, mu = 0)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns quantiles for for inverse Gaussian DF

    # FUNCTION:

    # Quantiles:
    # qgh(p = p, alpha = alpha, beta = beta, delta = delta, mu = mu,
    #   lambda = -0.5)

    # Checks:
    if (alpha <= 0) stop("alpha must be greater than zero")
    if (delta <= 0) stop("delta must be greater than zero")
    if (abs(beta) >= alpha) stop("abs value of beta must be less than alpha")

    # Internal Function:
    .froot <- function(x, alpha, beta, delta, mu, p)
        pnig(q = x, alpha = alpha, beta = beta, delta = delta, mu = mu) - p

    # Quantiles:
    ans <- p
    for(i in seq_along(p)) {
        lower = -1
        upper = +1
        counter = 0
        v <- NA
        while(is.na(v) && counter < 1000) {
            v <- .unirootNA(f = .froot, interval = c(lower, upper),
                            alpha = alpha, beta = beta, delta = delta,
                            mu = mu, p = p[i])
            counter = counter + 1
            lower = lower - 2^counter
            upper = upper + 2^counter
        }
        ans[i] <- v
    }

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


rnig <-
function(n, alpha = 1, beta = 0, delta = 1, mu = 0)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Return normal inverse Gaussian distributed random variates

    # Arguments:
    #   n - number of deviates to be generated
    #   alpha, beta - Shape Parameter, |beta| <= alpha
    #   delta  - Scale Parameter, 0 <= delta
    #   mu - Location Parameter

    # FUNCTION:

    # Settings:
    gamma = sqrt(alpha*alpha - beta*beta)

    # GAMMA:
    if (gamma == 0) {
        # GAMMA = 0:
        V = rnorm(n)^2
        Z = delta*delta / V
        X = sqrt(Z)*rnorm(n)
    } else {
        # GAMMA > 0:
        U = runif(n)
        V = rnorm(n)^2
        # FIXED ...
        z1 <- function(v, delta, gamma) {
            delta/gamma + v/(2*gamma^2) - sqrt(v*delta/(gamma^3) +
            (v/(2*gamma^2))^2 )
        }
        z2 <- function(v, delta, gamma) {
            (delta/gamma)^2 / z1(v = v, delta = delta, gamma = gamma)
        }
        pz1 <- function(v, delta, gamma) {
            delta / (delta + gamma * z1(v = v, delta = delta, gamma = gamma) )
        }
        s = (1-sign(U-pz1(v = V, delta = delta, gamma = gamma)))/2
        Z = z1(v = V, delta = delta, gamma = gamma)*s + z2(v = V, delta =
            delta, gamma = gamma)*(1-s)
        X = mu + beta*Z + sqrt(Z)*rnorm(n)
    }

    # Attributes:
    attr(X, "control") = c(dist = "nig", alpha = alpha, beta = beta,
        delta = delta, mu = mu)


    # Return Value:
    X
}


################################################################################


.qnigC <-
function(p, alpha = 1, beta = 0, delta = 1, mu = 0)
{
    # Description:
    #   Returns quantiles for for inverse Gaussian DF

    # Example:
    #   p = runif(10); .qnigC(p)

    # FUNCTION:

    # Checks:
    if(alpha <= 0) stop("Invalid parameters: alpha <= 0.0\n")
    if(alpha^2 <= beta^2) stop("Invalid parameters: alpha^2 <= beta^2.0\n")
    if(delta <= 0) stop("Invalid parameters: delta <= 0.0\n")
    if((sum(is.na(p)) > 0))
        stop("Invalid probabilities:\n",p,"\n")
    else
        if(sum(p < 0)+sum(p > 1) > 0) stop("Invalid probabilities:\n",p,"\n")

    # Evaluate NIG cdf by calling C function
    n <- length(p)
    q <- rep(0, n)
    retValues <- .C("qNIG",
        as.double(p),
        as.double(mu),
        as.double(delta),
        as.double(alpha),
        as.double(beta),
        as.integer(n),
        as.double(q),
        PACKAGE = "fBasics")
    quantiles <- retValues[[7]]
    quantiles[quantiles <= -1.78e+308] <- -Inf
    quantiles[quantiles >= 1.78e+308] <- Inf

    # Return Value:
    quantiles
}


# ------------------------------------------------------------------------------


.pnigC <-
function(q, alpha = 1, beta = 0, delta = 1, mu = 0)
{
    # Description:
    #   Returns probabilities for for inverse Gaussian DF

    # IMPORTANT NOTE:
    #   DW: C program fails, remains to check

    # Example:
    #   q = sort(rnorm(10)); .pnigC(q) # FAILS

    # FUNCTION:

    # Checks:
    if(alpha <= 0) stop("Invalid parameters: alpha <= 0.0\n")
    if(alpha^2 <= beta^2) stop("Invalid parameters: alpha^2 <= beta^2.0\n")
    if(delta <= 0) stop("Invalid parameters: delta <= 0.0\n")
    if(sum(is.na(q)) > 0) stop("Invalid quantiles:\n", q)

    # Evaluate NIG cdf by calling C function
    n <- length(q)
    p <- rep(0, n)
    retValues <- .C("pNIG",
        as.double(q),
        as.double(mu),
        as.double(delta),
        as.double(alpha),
        as.double(beta),
        as.integer(n),
        as.double(p),
        PACKAGE = "fBasics")
    probs <- retValues[[7]]

    # Return Value:
    probs
}


################################################################################

