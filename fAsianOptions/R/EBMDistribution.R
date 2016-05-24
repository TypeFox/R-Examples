
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
#   1999 - 2004, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file 


################################################################################
# FUNCTION:           EBM DENSITY APPROXIMATIONS:
#  dlognorm            log-Normal density an derivatives
#  plognorm            log-Normal, synonyme for plnorm
#  dgam                Gamma density, synonyme for dgamma
#  pgam                Gamma probability, synonyme for pgamma
#  drgam               Reciprocal-Gamma density
#  prgam               Reciprocal-Gamma probability
#  djohnson            Johnson Type I density
#  pjohnson            Johnson Type I probability
# FUNCTION :          MOMENTS FOR EBM DENSITY APPROXIMATIONS:
#  mnorm               Moments of Normal density
#  mlognorm            Moments of log-Normal density
#  mrgam               Moments of reciprocal-Gamma density
#  mjohnson            Moments of the Johnson Type-I density
#  masian              Moments of Asian Option density
#  .DufresneMoments     Internal Function used by masian()
# FUNCTION:           NUMERICAL DERIVATIVES:
#  derivative          First and second numerical derivative
# FUNCTION:           ASIAN DENSITY:
#  d2EBM               Double Integrated EBM density
#  .thetaEBM            Internal Function used to compute *2EBM()
#  .psiEBM              Internal Function used to compute *2EBM()
#  dEBM                Exponential Brownian motion density
#  pEBM                Exponential Brownian motion probability              
#  .gxuEBM              Internal Function used to compute *EBM()
#  .gxtEBM              Internal Function used to compute *EBM()
#  .gxtuEBM             Internal Function used to compute *EBM()
#  dasymEBM            Exponential Brownian motion asymptotic density
################################################################################


################################################################################
# EBM DENSITY APPROXIMATIONS:

dlognorm = 
function(x, meanlog = 0, sdlog = 1, deriv = c(0, 1, 2))
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates the log-Normal density or its first or
    #   second derivative. 
    
    # Arguments:
    
    # Details:
    #   Uses the function dlnorm().
    
    # See also:
    #   dlnorm(x, meanlog = 0, sdlog = 1, log = FALSE)
    #   plnorm(q, meanlog = 0, sdlog = 1, lower.tail = TRUE, log.p = FALSE)
    #   qlnorm(p, meanlog = 0, sdlog = 1, lower.tail = TRUE, log.p = FALSE)
    #   rlnorm(n, meanlog = 0, sdlog = 1)
    
    # FUNCTION:
    
    # Settings:
    deriv = deriv[1]
    
    # Function:
    result = dlnorm(x, meanlog = meanlog, sdlog = sdlog) 
    
    # First derivative, if desired:
    if (deriv == 1) {
        h1 = -(1/x + (log(x)-meanlog)/(sdlog^2*x))
        result = result * h1 
    }
        
    # Second derivative, if desired:    
    if (deriv == 2) {
        h1 = -(1/x + (log(x)-meanlog)/(sdlog^2*x))
        h2 = -(-1/x^2 + (-1/x^2)*(log(x)-meanlog)/sdlog^2 + 1/(sdlog^2*x^2))
        result = result * (h1^2 + h2) 
    }
        
    # Return Value:
    result
}


# ------------------------------------------------------------------------------


plognorm = 
function(q, meanlog = 0, sdlog = 1)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates the log-Normal probability. 
    
    # Details:
    #   Uses the function plnorm().
    
    # See also:
    #   dlnorm(x, meanlog = 0, sdlog = 1, log = FALSE)
    #   plnorm(q, meanlog = 0, sdlog = 1, lower.tail = TRUE, log.p = FALSE)
    #   qlnorm(p, meanlog = 0, sdlog = 1, lower.tail = TRUE, log.p = FALSE)
    #   rlnorm(n, meanlog = 0, sdlog = 1)
    
    # FUNCTION:
    
    # Resul:
    result = plnorm(q = q, meanlog = meanlog, sdlog = sdlog) 

    # Return Value:
    result
}
    

# ******************************************************************************


dgam = 
function(x, alpha, beta)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates the Gamma density.

    # Details:
    #   The function is a synonym to "dgamma".
    
    # See also:
    #   dgamma(x, shape, rate=1, scale=1/rate, log = FALSE)
    #   pgamma(q, shape, rate=1, scale=1/rate, lower.tail = TRUE, log = FALSE)
    #   qgamma(p, shape, rate=1, scale=1/rate, lower.tail = TRUE, log = FALSE)
    #   rgamma(n, shape, rate=1, scale=1/rate)

    # FUNCTION:
    
    # Return Value:
    dgamma(x = x, shape = alpha, scale = beta)
}


# ------------------------------------------------------------------------------


pgam = 
function(q, alpha, beta, lower.tail = TRUE)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates the Gamma probability.

    # Details:
    #   The function is a synonym to "pgamma".
    
    # See also:
    #   dgamma(x, shape, rate=1, scale=1/rate, log = FALSE)
    #   pgamma(q, shape, rate=1, scale=1/rate, lower.tail = TRUE, log.p = FALSE)
    #   qgamma(p, shape, rate=1, scale=1/rate, lower.tail = TRUE, log.p = FALSE)
    #   rgamma(n, shape, rate=1, scale=1/rate)

    # FUNCTION:
    
    # Return Value:
    pgamma(q = q, shape = alpha, scale = beta, lower.tail = lower.tail)
}


# ------------------------------------------------------------------------------


drgam = 
function(x, alpha = 1, beta = 1, deriv = c(0, 1, 2))
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates the reciprocal-Gamma density.

    # See also:
    #   dgamma(x, shape, rate=1, scale=1/rate, log = FALSE)
    #   pgamma(q, shape, rate=1, scale=1/rate, lower.tail = TRUE, log.p = FALSE)
    #   qgamma(p, shape, rate=1, scale=1/rate, lower.tail = TRUE, log.p = FALSE)
    #   rgamma(n, shape, rate=1, scale=1/rate)

    # FUNCTION:
    
    # Function Value:
    deriv = deriv[1]
    gr = dgamma(x = 1/x, shape = alpha, scale = beta) / (x^2)
    result = gr
    
    # First Derivative:
    if (deriv == 1) {
        h = function(x, alpha, beta) { 
            -(alpha+1)/x + 1/(beta*x^2) 
        }
        gr1 = gr*h(x, alpha, beta)
        result = gr1 
    }
    
    # Second Derivative:
    if (deriv == 2) {
        h = function(x, alpha, beta) { 
            -(alpha+1)/x + 1/(beta*x^2)
        }
        h1 = function(x, alpha, beta) { 
            +(alpha+1)/x^2 - 2/(beta*x^3) 
        }
        gr2 = gr*(h(x, alpha, beta)^2 + h1(x, alpha, beta))
        result = gr2 
    }
    
    # Return Value:
    result
}


# ------------------------------------------------------------------------------


prgam = 
function(q, alpha = 1, beta = 1, lower.tail = TRUE)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates the reciprocal-Gamma probability
    
    # See also:
    #   dgamma(x, shape, rate=1, scale=1/rate, log = FALSE)
    #   pgamma(q, shape, rate=1, scale=1/rate, lower.tail = TRUE, log.p = FALSE)
    #   qgamma(p, shape, rate=1, scale=1/rate, lower.tail = TRUE, log.p = FALSE)
    #   rgamma(n, shape, rate=1, scale=1/rate)

    # FUNCTION:
    
    # Return Value:
    1 - pgamma(q = 1/q, shape = alpha, scale = beta, lower.tail = lower.tail)
}


# ------------------------------------------------------------------------------


djohnson = 
function(x, a = 0, b = 1, c = 0, d = 1, deriv = c(0, 1, 2))
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates the Johnson Type-I density.
    
    # FUNCTION:
    
    # Function Value:
    deriv = deriv[1]
    z = a + b * log( (x-c)/d )
    z1 = b / (x-c)
    phi = dnorm(z, mean = 0, sd = 1)
    johnson = phi * z1
    result = johnson
    
    # First Derivative:
    if (deriv == 1) {
        z2 = -b / (x-c)^2
        johnson1 = phi * ( z2 - z*z1^2 )
        result = johnson1 
    }
    
    # Second Derivative:
    if (deriv == 2) {
        z2 = -b / (x-c)^2
        z3 = 2 * b / (x-c)^3
        johnson2 = phi * ( - z*z1*z2 + z^2*z1^3 + z3 -z1^3 - 2*z*z1*z2 )
        result = johnson2 
    }
    
    # Return Value:
    result
}


# ------------------------------------------------------------------------------


pjohnson = 
function(q, a = 0, b = 1, c = 0, d = 1)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates the Johnson Type-I probability.
    
    # FUNCTION:
    
    # Type I:
    z = a + b * log( (q-c) / d )
    
    # Return Value:
    pnorm(q = z, mean = 0, sd = 1)
}


################################################################################
# MOMENTS FOR EBM DENSITY APPROXIMATIONS:


mnorm = 
function(mean = 0, sd = 1)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes the moments for the Normal distribution.
    
    # FUNCTION:
    
    # Raw Moments:
    M = c( 
        mean, mean^2+sd^2, mean*(mean^2+3*sd*2), mean^4+6*mean^2*sd^2+3*sd^4 )
    
    # Centered Moments:
    m = M
    m[2] = M[2] - 2*M[1]*M[1] +   M[1]^2 
    m[3] = M[3] - 3*M[2]*M[1] + 3*M[1]*M[1]^2 -   M[1]^3 
    m[4] = M[4] - 4*M[3]*M[1] + 6*M[2]*M[1]^2 - 4*M[1]*M[1]^3 + M[1]^4 
    
    # Fischer Parameters - Skewness and Kurtosis:
    f = c(NA, NA)
    f[1] = m[3] / m[2]^(3/2)
    f[2] = m[4] / m[2]^2 - 3
    
    # Return Value:
    list(rawMoments = M, centralMoments = m, fisher = f)
}


# ------------------------------------------------------------------------------


mlognorm = 
function(meanlog = 0, sdlog = 1)
{   # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Computes the moments for the Log-Normal distribution.
    
    # FUNCTION:
    
    # Raw Moments:
    n = 1:4
    M = exp ( n * meanlog + n^2 * sdlog^2/2 )
    
    # Centered Moments:
    m = M
    m[2] = M[2] - 2*M[1]*M[1] +   M[1]^2 
    m[3] = M[3] - 3*M[2]*M[1] + 3*M[1]*M[1]^2 -   M[1]^3 
    m[4] = M[4] - 4*M[3]*M[1] + 6*M[2]*M[1]^2 - 4*M[1]*M[1]^3 + M[1]^4 
    
    # Fischer Parameters - Skewness and Kurtosis:
    f = c(NA, NA)
    f[1] = m[3] / m[2]^(3/2)
    f[2] = m[4] / m[2]^2 - 3
    
    # Return Value:
    list(rawMoments = M, centralMoments = m, fisher = f)
}


# ------------------------------------------------------------------------------


mrgam = 
function(alpha = 1/2, beta = 1)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes the moments for the Reciprocal-Gamma distribution.
    
    # FUNCTION:
    
    # Raw Moments:
    M = rep(0, times = 4)
    M[1] =   1  / (beta*(alpha - 1))
    M[2] = M[1] / (beta*(alpha - 2))
    M[3] = M[2] / (beta*(alpha - 3))
    M[4] = M[3] / (beta*(alpha - 4))
    
    # Centered Moments:
    m = M
    m[2] = M[2] - 2*M[1]*M[1] +   M[1]^2 
    m[3] = M[3] - 3*M[2]*M[1] + 3*M[1]*M[1]^2 -   M[1]^3 
    m[4] = M[4] - 4*M[3]*M[1] + 6*M[2]*M[1]^2 - 4*M[1]*M[1]^3 + M[1]^4 
    
    # Fischer Parameters - Skewness and Kurtosis:
    f = c(NA, NA)
    f[1] = m[3] / m[2]^(3/2)
    f[2] = m[4] / m[2]^2 - 3
    
    # Return Value:
    list(rawMoments = M, centralMoments = m, fisher = f)
}


# ------------------------------------------------------------------------------


mjohnson = 
function(a, b, c, d)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes the moments for the Johnson Type-I distribution
    
    # FUNCTION:
    
    # Raw Moments:
    M = c(NA, NA, NA, NA)
    
    # Centered Moments:
    m = M
    m[2] = M[2] - 2*M[1]*M[1] +   M[1]^2 
    m[3] = M[3] - 3*M[2]*M[1] + 3*M[1]*M[1]^2 -   M[1]^3 
    m[4] = M[4] - 4*M[3]*M[1] + 6*M[2]*M[1]^2 - 4*M[1]*M[1]^3 + M[1]^4 
    
    # Fischer Parameters - Skewness and Kurtosis:
    f = c(NA, NA)
    f[1] = m[3] / m[2]^(3/2)
    f[2] = m[4] / m[2]^2 - 3
    
    # Return Value:
    list(rawMoments = M, centralMoments = m, fisher = f)
}


# ------------------------------------------------------------------------------


masian = 
function(Time = 1, r = 0.045, sigma = 0.30)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes the moments for the Asian-Option distribution
    
    # FUNCTION:
    
    # Raw Moments:
    M = .DufresneMoments(M = 4, Time = Time, r = r, sigma = sigma)

    # Centered Moments:
    m = M
    m[2] = M[2] - 2*M[1]*M[1] +   M[1]^2 
    m[3] = M[3] - 3*M[2]*M[1] + 3*M[1]*M[1]^2 -   M[1]^3 
    m[4] = M[4] - 4*M[3]*M[1] + 6*M[2]*M[1]^2 - 4*M[1]*M[1]^3 + M[1]^4 
    
    # Fischer Parameters - Skewness and Kurtosis:
    f = c(NA, NA)
    f[1] = m[3] / m[2]^(3/2)
    f[2] = m[4] / m[2]^2 - 3
    
    # Return Value:
    list(rawMoments = M, centralMoments = m, fisher = f)
}


# ------------------------------------------------------------------------------


.DufresneMoments = 
function (M = 4, Time = 1, r = 0.045, sigma = 0.30) 
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes raw moments for the Asian-Option distribution.
    
    # Note:
    #   Called by function masian()
    
    # FUNCTION:
    
    # Internal Function:
    moments = 
    function (M, tau, nu) { 
        d = function(j, n, beta) {
            d = 2^n
            for (i in 0:n) if (i != j) d = d  / ( (beta+j)^2 - (beta+i)^2 )
            d 
        }     
        moments = rep(0, length = M)
        for (n in 1:M) {    
            moments[n] = 0
            for (j in 0:n) moments[n] = moments[n] + 
                d(j, n, nu/2)*exp(2*(j^2+j*nu)*tau)
            moments[n] = prod(1:n) * moments[n] / (2^(2*n)) 
        }
        moments 
    }
    
    # Compute:
    tau = sigma^2*Time/4
    nu = 2*r/sigma^2-1
    ans = (4/sigma^2)^(1:M) * moments(M, tau, nu)
    
    # Return Value:
    ans 
}
        

################################################################################


derivative = 
function(x, y, deriv = c(1, 2))
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates numerically the first or second derivative
    #   of the functuin y(x) by finite differences.
    
    # Arguments:
    #   x - a numeric vector of values
    #   y - a numeric vectror of function values y(x)
    #   deriv - the degree of differentiation, either 1 or 2.
    
    # FUNCTION:
    
    # Stop in the case of wrong argument deriv:
    dseriv = deriv[1]
    if (deriv < 1 || deriv > 2) stop("argument error")
    
    # Function to calculate the next derivative by differences:
    "calcderiv" = function(x,y) {
        list(x = x[2:length(x)]-diff(x)/2, y = diff(y)/diff(x))}
    
    # First Numerical Derivative:
    result = calcderiv(x, y)
    
    # Second Numerical Derivative, if desired:
    if (deriv == 2) result = calcderiv(result$x, result$y)
    
    # Return Value:
    list(x = result$x, y = result$y)
}


################################################################################


d2EBM =
function(u, t = 1) 
{   # A function written by Diethelm Wuertz

    # Description:
    #   Calculate the density integral "f_A_t(u)" given by  
    #   equation 4.36 in: R. Gould, "The Distribution of the 
    #   Integral of Exponential Brownian Motion".
    
    # Arguments:
    #   t - numeric value
    #   u - numeric value
    
    # FUNCTION:
    
    # Function to be integrated:
    f = function(x, tt, uu) {
        fx = rep(0, length=length(x))
        for (i in 1:length(x) )
            fx[i] = (1/uu) * exp(-(1+exp(2*x[i]))/(2*uu)) *  
                .thetaEBM(r=exp(x[i])/uu, u=tt) 
        fx }
    
    # Integrate:
    result = rep(0, length = length(u))
    for (i in 1:length(u)) {
    result[i] = integrate(f, lower = -16, upper = 4, tt = t, uu = u[i], 
        subdivisions = 100, rel.tol=.Machine$double.eps^0.25, 
        abs.tol=.Machine$double.eps^0.25)$value }
    
    # Return Value:
    result
}


# ------------------------------------------------------------------------------


.thetaEBM =
function(r, u) 
{   # A function written by Diethelm Wuertz

    # Description:
    #   Calculate the integral "\theta_r(u)" given by equations 
    #   2.22 and 2.23 in: R. Gould, "The Distribution of the 
    #   Integral of Exponential Brownian Motion".
    
    # Arguments:
    #   r - vector of numeric values
    #   u - numeric value
    
    # FUNCTION:
    
    # Function to be integrated:  
    f = function(x, rr, uu) {
        fx = rep(0, length = length(x))
        for (i in 1:length(x) ) 
            fx[i] = exp(-x[i]^2/(2*uu)) * exp(-rr*cosh(x[i])) * 
                sinh(x[i]) * sin(pi*x[i]/uu) 
        fx }
    
    # Loop over r-Vector:
    result = rep(0, length=length(r))
    for ( i in 1: length(r) ) {
        result[i] = integrate(f, lower = 0, upper = 30, rr = r[i], uu = u, 
            subdivisions = 100, rel.tol = .Machine$double.eps^0.25, 
            abs.tol = .Machine$double.eps^0.25)$value 
        result[i] = result[i] * (r/sqrt((2*u*pi^3))) * exp(pi^2/(2*u)) }
    
    # Return Value:
    result
}


# ------------------------------------------------------------------------------


.psiEBM = 
function(r, u)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Calculate the integral "\psi_r(u)" given by equations 
    #   2.22 and 2.23 in: R. Gould, "The Distribution of the 
    #   Integral of Exponential Brownian Motion".
    
    # Arguments:
    #   r - vector of numeric values
    #   u - numeric value
    
    # FUNCTION:
    
    # Calculate psi() from theta():
    result = sqrt(2*u*pi^3) * exp(-pi^2/(2*u)) * .thetaEBM(r, u)
    
    # Return Value:
    result
    
}


# ------------------------------------------------------------------------------


dEBM =
function(u, t = 1)
{   # A function written by Diethelm Wuertz

    # Arguments;
    #   t - a numeric value
    #   u - a vector of numeric values
    
    # FUNCTION:
    
    # Calculate Density:
    result = rep(0, times = length(u))
    for (i in 1:length(u) ) {
        result[i] = integrate(.gxtuEBM, lower = 0, upper = 100, 
            t = t, u = u[i])$value 
    }
        
    # Return Value:
    result
}


# ------------------------------------------------------------------------------


pEBM = 
function(u, t = 1)
{   # A function written by Diethelm Wuertz

    # Arguments;
    #   t - a numeric value
    #   u - a vector of numeric value
    
    # FUNCTION:
    
    # Calculate Probability:
    result = rep(0, times = length(u))
    result[1] = integrate(dEBM, lower = 0, upper = u[1], t = t)$value
    if (length(u) > 1) {
        for (i in 2:length(u) ) {
            result[i] = result[i-1] + integrate(
                dEBM, lower = u[i-1], upper = u[i], t = t)$value 
        } 
    }
        
    # Return Value:
    result
}


# ------------------------------------------------------------------------------


.gxuEBM =
function(x, u)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Interchange the Integrals - and first integrate:
    #   1/u^2 * sinh(x) * exp(-(1/(2*u))*(1+exp(2*x))) * exp(x) * 
    #       exp(-exp(x)*cosh(y)/u)
    
    # FUNCTION:
    
    # Compute g(x, u): 
    fx = rep(0, length = length(x)) 
    if (u > 0) {
        for ( i in 1:length(x) ) {
            su = (u)^(-3/2) * sinh(x[i])
            cx = cosh(x[i])/sqrt(u)
            sx = sinh(x[i])/sqrt(u)
            Asymptotics = exp(x[i]) / sqrt(u) / 2
            if (Asymptotics <= 33.0) {
                fx[i] = su * pnorm(-cx) / dnorm(sx) }
            else {
                fx[i] = su * exp(-1/(2*u)) * (1-1/cx^2+3/cx^4) / cx  } } } 
                
    # Return Value:
    fx
}


# ------------------------------------------------------------------------------


.gxtEBM = 
function(x, t)
{   # A function written by Diethelm Wuertz

    # Description:
    
    # FUNCTION:
    
    # Compute g(x, t):
    fx = exp(pi^2/(2*t)) / sqrt(2*t*pi^3) * exp(-x^2/(2*t))

    # Return Value:
    fx
}


# ------------------------------------------------------------------------------


.gxtuEBM =  
function(x, t, u)
{   # A function written by Diethelm Wuertz

    # FUNCTION:
    
    # Result:
    fx = .gxtEBM(x = x, t = t) * .gxuEBM(x = x, u = u) * sin(pi*x/t) 

    # Return Value:
    fx
}


# ------------------------------------------------------------------------------


dasymEBM = 
function(u, t = 1)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Calculates the asymptotic behavior of the density
    #   function f of the exponential Brownian maotion
    
    # FUNCTION:
    
    # Asymptotic Density:
    alpha = log ( 8*u*exp(-2*t) )
    beta = exp (  -((log(alpha/(4*t)))^2)/(8*t)  )
    f = sqrt(t) * exp(t/2) * exp(-alpha^2/(8*t)) * beta
    
    # Take care of gamma function ...
    warn = options()$warn
    options(warn = -1)  
    # f = f / (u * sqrt(u) * alpha * gamma(alpha/(4 * t)))
    f = f / (u * sqrt(u) * (4*t) * gamma(alpha/(4 * t)+1))
    f[is.na(f)] = 0
    options(warn = warn)  
    
    # Return Value:
    f
}


################################################################################

