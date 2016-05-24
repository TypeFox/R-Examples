
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
# FUNCTION:             
#  .rghyp   
#  .rgigjd    
#  .rgigjd1
#  .dhyp1
#  .dhyp2
#  .dhyp3
#  .dhyp4
#  .phyp1
#  .phyp2
#  .phyp3
#  .phyp4
#  .qhyp1
#  .qhyp2
#  .qhyp3
#  .qhyp4
#  .rhyp1
#  .rhyperb
#  .hyperb.change.pars
#  .rhyp2
#  .rhyp3
#  .rhyp4
################################################################################


# This is code borrowed from 
#   David Scott's website.
#   The functions are also part of previous versions of the R contributed
#   package Hyperbolic Dist.


# Rmetrics:
#   Note that these functions are not packaged and available on Debian 
#   as of 2007-06-23. 
#   To run these functions under Debian/Rmetrics we have them    
#   implemented here as a builtin.


################################################################################


.rghyp <-
function(n, theta)
{
    # A function implemented by Diethelm Wuertz

    # Author:
    #   Original Version by David Scott

    # FUNCTION:

    # Settings:
    lambda = theta[1]
    alpha = theta[2]
    beta = theta[3]
    delta = theta[4]
    mu = theta[5]
    chi = delta^2
    psi = alpha^2 - beta^2

    # Ckecks:
    if (alpha <= 0) stop("alpha must be greater than zero")
    if (delta <= 0) stop("delta must be greater than zero")
    if (abs(beta) >= alpha) stop("abs value of beta must be less than alpha")

    # Random Numbers:
    if (lambda == 1){
        X = .rgigjd1(n, c(lambda, chi, psi))
    } else{
        X = .rgigjd(n, c(lambda, chi, psi))
    }

    # Result:
    sigma = sqrt(X)
    Z = rnorm(n)
    Y = mu + beta*sigma^2 + sigma*Z

    # Return Value:
    Y
}


# ------------------------------------------------------------------------------


.rgigjd <-
function(n, theta)
{
    # A function implemented by Diethelm Wuertz

    # Author:
    #   Original Version by David Scott

    # FUNCTION:

    # Settings:
    lambda = theta[1]
    chi = theta[2]
    psi = theta[3]

    # Checks:
    if (chi < 0) stop("chi can not be negative")
    if (psi < 0) stop("psi can not be negative")
    if ((lambda >= 0)&(psi==0)) stop("When lambda >= 0, psi must be > 0")
    if ((lambda <= 0)&(chi==0)) stop("When lambda <= 0, chi must be > 0")
    if (chi == 0) stop("chi = 0, use rgamma")
    if (psi == 0) stop("algorithm only valid for psi > 0")

    alpha = sqrt(psi/chi)
    beta = sqrt(psi*chi)

    m = (lambda-1+sqrt((lambda-1)^2+beta^2))/beta

    g = function(y){
        0.5*beta*y^3 - y^2*(0.5*beta*m+lambda+1) +
            y*((lambda-1)*m-0.5*beta) + 0.5*beta*m
    }

    upper = m
    while (g(upper) <= 0) upper = 2*upper
    yM = uniroot(g, interval=c(0,m))$root
    yP = uniroot(g, interval=c(m,upper))$root

    a = (yP-m)*(yP/m)^(0.5*(lambda-1))*exp(-0.25*beta*(yP+1/yP-m-1/m))
    b = (yM-m)*(yM/m)^(0.5*(lambda-1))*exp(-0.25*beta*(yM+1/yM-m-1/m))
    c = -0.25*beta*(m+1/m) + 0.5*(lambda-1)*log(m)

    output = numeric(n)

    for(i in 1:n){
        need.value = TRUE
        while(need.value==TRUE){
            R1 = runif (1)
            R2 = runif (1)
            Y = m + a*R2/R1 + b*(1-R2)/R1
            if (Y>0){
                if (-log(R1)>=-0.5*(lambda-1)*log(Y)+0.25*beta*(Y+1/Y)+c){
                    need.value = FALSE
                }
            }
        }
        output[i] = Y
    }

    # Return Value:
    output/alpha
}


# ------------------------------------------------------------------------------


.rgigjd1 <-
function(n, theta)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Modified version of rgigjd to generate random observations
    #   from a generalised inverse Gaussian distribution in the
    #   special case where lambda = 1.

    # Author:
    #   Original Version by David Scott

    # FUNCTION:

    if (length(theta) == 2) theta = c(1, theta)

    # Settings:
    lambda = 1
    chi = theta[2]
    psi = theta[3]

    # Checks:
    if (chi < 0) stop("chi can not be negative")
    if (psi < 0) stop("psi can not be negative")
    if (chi == 0) stop("chi = 0, use rgamma")
    if (psi == 0) stop("When lambda >= 0, psi must be > 0")

    alpha = sqrt(psi/chi)
    beta = sqrt(psi*chi)
    m = abs(beta)/beta
    g = function(y){
        0.5*beta*y^3 - y^2*(0.5*beta*m+lambda+1) +
            y*(-0.5*beta) + 0.5*beta*m
    }

    upper = m
    while (g(upper)<=0) upper = 2*upper
    yM = uniroot(g,interval=c(0,m))$root
    yP = uniroot(g,interval=c(m,upper))$root

    a = (yP-m)*exp(-0.25*beta*(yP+1/yP-m-1/m))
    b = (yM-m)*exp(-0.25*beta*(yM+1/yM-m-1/m))
    c = -0.25*beta*(m+1/m)

    output = numeric(n)

    for(i in 1:n){
        need.value = TRUE
        while(need.value==TRUE){
            R1 = runif (1)
            R2 = runif (1)
            Y = m + a*R2/R1 + b*(1-R2)/R1
            if (Y>0){
                if (-log(R1)>=0.25*beta*(Y+1/Y)+c){
                    need.value = FALSE
                }
            }
        }
        output[i] = Y
    }

    # Return Value:
    output/alpha
}


# ------------------------------------------------------------------------------


.dhyp1 <-
function(x, alpha = 1, beta = 0, delta = 1, mu = 0)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns Hyperbolic Density Function PDF

    # Arguments:
    #   alpha, beta - Shape Parameter, |beta| <= alpha
    #   delta  - Scale Parameter, 0 <= delta
    #   mu - Location Parameter

    # FUNCTION:

    # Density:
    efun = exp( -alpha*sqrt(delta^2 + (x-mu)^2) + beta*(x-mu) )
    sqr = sqrt(alpha^2-beta^2)
    bK1 = besselK(delta*sqr, nu = 1)
    prefac = sqr / ( 2 * alpha * delta * bK1)
    ans = prefac * efun

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.dhyp2 <-
function(x, zeta = 1, rho = 0, delta = 1, mu = 0)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns Hyperbolic density in the 2nd parameterization

    # FUNCTION:

    # Arguments:
    #   alpha, beta - Shape Parameter, |beta| <= alpha
    #   delta  - Scale Parameter, 0 <= delta
    #   mu - Location Parameter

    # Parameter Change:
    alpha = zeta / ( delta * sqrt(1 - rho*rho) )
    beta = alpha * rho

    # Return Value:
    ans = dhyp(x, alpha, beta, delta, mu)
    ans
}


# ------------------------------------------------------------------------------


.dhyp3 <-
function(x, xi = 1/sqrt(2), chi = 0, delta = 1, mu = 0)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns Hyperbolic density in the 2nd parameterization

    # Arguments:
    #   alpha, beta - Shape Parameter, |beta| <= alpha
    #   delta  - Scale Parameter, 0 <= delta
    #   mu - Location Parameter

    # FUNCTION:

    # Parameter Change:
    rho = chi / xi
    zeta = 1/xi^2 - 1
    alpha = zeta / ( delta * sqrt(1 - rho*rho) )
    beta = alpha * rho

    # Return Value:
    ans = dhyp(x, alpha, beta, delta, mu)
    ans
}


# ------------------------------------------------------------------------------


.dhyp4 <-
function(x, a.bar = 1, b.bar = 0, delta = 1, mu = 0)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns Hyperbolic density in the 2nd parameterization

    # Arguments:
    #   alpha, beta - Shape Parameter, |beta| <= alpha
    #   delta  - Scale Parameter, 0 <= delta
    #   mu - Location Parameter

    # FUNCTION:

    # Parameter Change:
    alpha = a.bar / delta
    beta = b.bar / delta

    # Return Value:
    ans = dhyp(x, alpha, beta, delta, mu)
    ans
}


# ------------------------------------------------------------------------------


.phyp1 <-
function(q, alpha = 1, beta = 0, delta = 1, mu = 0, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Return cumulative probability of Hyperbolic PDF

    # Arguments:
    #   alpha, beta - Shape Parameter, |beta| <= alpha
    #   delta  - Scale Parameter, 0 <= delta
    #   mu - Location Parameter

    # FUNCTION:

    # Cumulative Probability:
    ans = NULL
    for (Q in q) {
        Integral = integrate(dhyp, -Inf, Q, stop.on.error = FALSE,
            alpha = alpha, beta = beta, delta = delta, mu = mu, ...)
        # Works in both, R and SPlus:
        ans = c(ans, as.numeric(unlist(Integral)[1]) )
    }

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.phyp2 <-
function(q, zeta = 1, rho = 0, delta = 1, mu = 0, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns cumulative probability in the 2nd parameterization

    # Arguments:
    #   zeta, rho - Shape Parameter, resulting |beta| <= alpha
    #   delta  - Scale Parameter, 0 <= delta
    #   mu - Location Parameter

    # FUNCTION:

    # Parameter Change:
    alpha = zeta / ( delta * sqrt(1 - rho*rho) )
    beta = alpha * rho

    # Return Value:
    ans = phyp(q, alpha, beta, delta, mu, ...)
    ans
}


# ------------------------------------------------------------------------------


.phyp3 <-
function(q, xi = 1/sqrt(2), chi = 0, delta = 1, mu = 0, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns cumulative probability in the 3rd parameterization

    # Arguments:
    #   xi, xhi - Shape Parameter, resulting |beta| <= alpha
    #   delta  - Scale Parameter, 0 <= delta
    #   mu - Location Parameter

    # FUNCTION:

    # Parameter Change:
    rho = chi / xi
    zeta = 1/xi^2 - 1
    alpha = zeta / ( delta * sqrt(1 - rho*rho) )
    beta = alpha * rho

    # Return Value:
    ans = phyp(q, alpha, beta, delta, mu, ...)
    ans
}


# ------------------------------------------------------------------------------


.phyp4 <-
function(q, a.bar = 1, b.bar = 0, delta = 1, mu = 0, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns cumulative probability in the 4th parameterization

    # Arguments:
    #   a.bar, b.bar - Shape Parameter, resulting |beta| <= alpha
    #   delta  - Scale Parameter, 0 <= delta
    #   mu - Location Parameter

    # FUNCTION:

    # Parameter Change:
    alpha = a.bar / delta
    beta = b.bar / delta

    # Return Value:
    ans = phyp(q, alpha, beta, delta, mu, ...)
    ans
}


# ------------------------------------------------------------------------------


.qhyp1 <-
function(p, alpha = 1, beta = 0, delta = 1, mu = 0, ...)
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

    # Internal Functions:
    .froot <-
    function(x, alpha, beta, delta, p)
        {
            phyp(q = x, alpha = alpha, beta = beta, delta = delta, mu = 0) - p
        }

    # Loop over all p's:
    result = NULL
    for (pp in p) {
        lower = -1
        upper = +1
        counter = 0
        iteration = NA
        while (is.na(iteration)) {
            iteration = .unirootNA(f = .froot, interval = c(lower, upper),
                alpha = alpha, beta = beta, delta = delta, p = pp)
            counter = counter + 1
            lower = lower-2^counter
            upper = upper+2^counter
        }
        result = c(result, iteration)
    }

    # Return Value:
    ans = result + mu
    ans
}


# ------------------------------------------------------------------------------


.qhyp2 <-
function(p, zeta = 1, rho = 0, delta = 1, mu = 0, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns quantiles of Hyperbolic PDF in the 2nd parameterization

    # Arguments:
    #   zeta, rho - Shape Parameter, resulting |beta| <= alpha
    #   delta  - Scale Parameter, 0 <= delta
    #   mu - Location Parameter

    # FUNCTION:

    # Parameter Change:
    alpha = zeta / ( delta * sqrt(1 - rho*rho) )
    beta = alpha * rho

    # Return Value:
    ans = qhyp(p, alpha, beta, delta, mu, ...)
    ans
}


# ------------------------------------------------------------------------------


.qhyp3 <-
function(p, xi = 1/sqrt(2), chi = 0, delta = 1, mu = 0, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns quantiles of Hyperbolic PDF in the 3rd parameterization

    # Arguments:
    #   zeta, chi - Shape Parameter, resulting |beta| <= alpha
    #   delta  - Scale Parameter, 0 <= delta
    #   mu - Location Parameter

    # FUNCTION:

    # Parameter Change:
    rho = chi / xi
    zeta = 1/xi^2 - 1
    alpha = zeta / ( delta * sqrt(1 - rho*rho) )
    beta = alpha * rho

    # Return Value:
    ans = qhyp(p, alpha, beta, delta, mu, ...)
    ans
}


# ------------------------------------------------------------------------------


.qhyp4 <-
function(p, a.bar = 1, b.bar = 0, delta = 1, mu = 0, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns quantiles of Hyperbolic PDF in the 4th parameterization

    # Arguments:
    #   a.bar, b.bar - Shape Parameter, resulting |beta| <= alpha
    #   delta  - Scale Parameter, 0 <= delta
    #   mu - Location Parameter

    # FUNCTION:

    # Parameter Change:
    alpha = NA
    beta = NA

    # Return Value:
    ans = qhyp(p, alpha, beta, delta, mu, ...)
    ans
}


# ------------------------------------------------------------------------------


.rhyp1 <-
function(n, alpha = 1, beta = 0, delta = 1, mu = 0)
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

    # Result - Use Standard Parameterization:
    Zeta = delta * sqrt(alpha^2 - beta^2)
    hyp.Pi = beta / sqrt(alpha^2 - beta^2)
    theta = c(hyp.Pi, Zeta, delta, mu)

    # Return Value:
    ans = .rhyperb(n = n, theta = theta)
    ans
}


# ------------------------------------------------------------------------------


.rhyperb <-
function(n, theta)
{
    # FUNCTION:

    # Internal Function:
    hyp.pi = theta[1]
    zeta = theta[2]
    delta = theta[3]
    mu = theta[4]
    alpha = as.numeric(.hyperb.change.pars(1, 2, theta))[1] * delta
    beta = as.numeric(.hyperb.change.pars(1, 2, theta))[2] * delta
    phi = as.numeric(.hyperb.change.pars(1, 3, theta))[1] * delta
    gamma = as.numeric(.hyperb.change.pars(1, 3, theta))[2] * delta
    theta.start = -sqrt(phi * gamma)
    t = -sqrt(gamma/phi)
    w = sqrt(phi/gamma)
    delta1 = exp(theta.start)/phi
    delta2 = (w - t) * exp(theta.start)
    delta3 = exp(-gamma * w)/gamma
    k = 1/(delta1 + delta2 + delta3)
    r = k * delta1
    v = 1 - k * delta3
    output = numeric(n)
    need.value = TRUE
    for (i in 1:n) {
        while (need.value == TRUE) {
            U = runif(1)
            E = rexp(1)
            if (U <= r) {
                x = 1/phi * log(phi * U/k)
                if (E >= alpha * (sqrt(1 + x^2) + x)) {
                  need.value = FALSE } }
            if ((U > r) & (U <= v)) {
                x = t - 1/phi + U * exp(-theta.start)/k
                if (E >= alpha * sqrt(1 + x^2) - beta * x + theta.start) {
                    need.value = FALSE
                }
            }
            if (U > v) {
                x = 1/gamma * log(k/gamma) - 1/gamma * log(1 - U)
                if (E >= alpha * (sqrt(1 + x^2) - x)) {
                    need.value = FALSE
                }
            }
        }
        output[i] = delta * x + mu
        need.value = TRUE
    }

    # Return Value:
    output
}


# ------------------------------------------------------------------------------


.hyperb.change.pars <-
function(from, to, theta)
{
    # FUNCTION:

    # Internal Function:
    delta <- theta[3]
    mu <- theta[4]
    hyperb.pi <- theta[1]
    zeta <- theta[2]
    if (from == 1 && to == 2) {
        alpha <- zeta * sqrt(1 + hyperb.pi^2)/delta
        beta <- zeta * hyperb.pi/delta
        output = c(alpha = alpha, beta = beta, delta = delta, mu = mu)
    }
    if (from == 1 && to == 3) {
        phi <- zeta/delta * (sqrt(1 + hyperb.pi^2) + hyperb.pi)
        gamma <- zeta/delta * (sqrt(1 + hyperb.pi^2) - hyperb.pi)
        output = c(phi = phi, gamma = gamma, delta = delta, mu = mu)
    }

    # Return Value:
    output
}


# ------------------------------------------------------------------------------


.rhyp2 <-
function(n, zeta = 1, rho = 0, delta = 1, mu = 0)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns random deviates of Hyperbolic PDF in the 2nd parameterization

    # Arguments:
    #   zeta, rho - Shape Parameter, resulting |beta| <= alpha
    #   delta  - Scale Parameter, 0 <= delta
    #   mu - Location Parameter

    # FUNCTION:

    # Parameter Change:
    alpha = zeta / ( delta * sqrt(1 - rho*rho) )
    beta = alpha * rho
    ans = rhyp(n, alpha, beta, delta, mu)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.rhyp3 <-
function(n, xi = 1/sqrt(2), chi = 0, delta = 1, mu = 0)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns random deviates of Hyperbolic PDF in the 3rd parameterization

    # Arguments:
    #   zeta, chi - Shape Parameter, resulting |beta| <= alpha
    #   delta  - Scale Parameter, 0 <= delta
    #   mu - Location Parameter

    # FUNCTION:

    # Parameter Change:
    rho = chi / xi
    zeta = 1/xi^2 - 1
    alpha = zeta / ( delta * sqrt(1 - rho*rho) )
    beta = alpha * rho
    ans = rhyp(n, alpha, beta, delta, mu)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.rhyp4 <-
function(n, a.bar = 1, b.bar = 0, delta  = 1, mu = 0)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns random deviates of Hyperbolic PDF in the 4th parameterization

    # Arguments:
    #   a.bar, b.bar - Shape Parameter, resulting |beta| <= alpha
    #   delta  - Scale Parameter, 0 <= delta
    #   mu - Location Parameter

    # FUNCTION:

    # Parameter Change:
    alpha = a.bar / delta
    beta = b.bar / delta
    ans = rhyp(n, alpha, beta, delta, mu)

    # Return Value:
    ans
}


################################################################################

