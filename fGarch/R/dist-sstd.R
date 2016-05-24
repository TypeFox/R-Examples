
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
# FUNCTION:             DESCRIPTION:
#  dsstd                 Density for the skewed Student-t Distribution
#  psstd                 Probability function for the skewed STD
#  qsstd                 Quantile function for the skewed STD
#  rsstd                 Random Number Generator for the skewed STD
# FUNCTION:             DESCRIPTION:
#  .dsstd                Internal, density for the skewed Student-t Distribution
#  .psstd                Internal, probability function for the skewed STD
#  .qsstd                Internal, quantile function for the skewed STD
#  .rsstd                Internal, random Number Generator for the skewed STD
################################################################################


dsstd <-
function(x, mean = 0, sd = 1, nu = 5, xi = 1.5, log = FALSE)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Compute the density function of the
    #   skewed Student-t distribution

    # FUNCTION:

    # Params:
    if (length(mean) == 4) {
        xi = mean[4]
        nu = mean[3]
        sd = mean[2]
        mean = mean[1]
    }  
    
    # Shift and Scale:
    result = .dsstd(x = (x-mean)/sd, nu = nu, xi = xi) / sd

    # Log:
    if(log) result = log(result)
    
    # Return Value:
    result
}


# ------------------------------------------------------------------------------


psstd <-
function(q, mean = 0, sd = 1, nu = 5, xi = 1.5)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Compute the distribution function of the
    #   skewed Student-t distribution

    # FUNCTION:

    # Shift and Scale:
    result = .psstd(q = (q-mean)/sd, nu = nu, xi = xi)

    # Return Value:
    result
}


# ------------------------------------------------------------------------------


qsstd <-
function(p, mean = 0, sd = 1, nu = 5, xi = 1.5)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Compute the quantile function of the
    #   skewed Student-t distribution

    # FUNCTION:

    # Shift and Scale:
    result = .qsstd(p = p, nu = nu, xi = xi) * sd + mean

    # Return Value:
    result
}


# ------------------------------------------------------------------------------


rsstd <-
function(n, mean = 0, sd = 1, nu = 5, xi = 1.5)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Generate random deviates from the
    #   skewed Student-t distribution

    # FUNCTION:

    # Shift and Scale:
    result = .rsstd(n = n, nu = nu, xi = xi) * sd + mean

    # Return Value:
    result
}


################################################################################


.dsstd <-
function(x, nu, xi)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Internal Function

    # FUNCTION:

    # For SPlus compatibility:
    if (!exists("beta"))
        beta <- function (a, b) exp( lgamma(a) + lgamma(b) -lgamma(a+b) )

    # Standardize:
    m1 = 2 * sqrt(nu-2) / (nu-1) / beta(1/2, nu/2)
    mu = m1*(xi-1/xi)
    sigma =  sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)
    z = x*sigma + mu

    # Compute:
    Xi = xi^sign(z)
    g = 2 / (xi + 1/xi)
    Density = g * dstd(x = z/Xi, nu = nu)

    # Return Value:
    Density * sigma
}


# ------------------------------------------------------------------------------


.psstd <-
function(q, nu, xi)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Internal Function

    # FUNCTION:

    # For SPlus compatibility:
    if (!exists("beta"))
        beta <- function (a, b) exp( lgamma(a) + lgamma(b) -lgamma(a+b) )

    # Standardize:
    m1 = 2 * sqrt(nu-2) / (nu-1) / beta(1/2, nu/2)
    mu = m1*(xi-1/xi)
    sigma =  sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)
    z = q*sigma + mu

    # Compute:
    Xi = xi^sign(z)
    g = 2 / (xi + 1/xi)
    Probability = Heaviside(z) - sign(z) * g * Xi * 
        pstd(q = -abs(z)/Xi, nu = nu)

    # Return Value:
    Probability
}


# ------------------------------------------------------------------------------


.qsstd <-
function(p, nu, xi)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Internal Function

    # FUNCTION:

    # For SPlus compatibility:
    if (!exists("beta"))
        beta <- function (a, b) exp( lgamma(a) + lgamma(b) -lgamma(a+b) )

    # Standardize:
    m1 = 2 * sqrt(nu-2) / (nu-1) / beta(1/2, nu/2)
    mu = m1*(xi-1/xi)
    sigma =  sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)

    # Compute:
    g = 2  / (xi + 1/xi)
    sig = sign(p-1/2)
    Xi = xi^sig
    p = (Heaviside(p-1/2)-sig*p) / (g*Xi)
    Quantile = (-sig*qstd(p = p, sd = Xi, nu = nu) - mu ) / sigma

    # Return Value:
    Quantile
}


# ------------------------------------------------------------------------------


.rsstd <-
function(n, nu, xi)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Internal Function

    # FUNCTION:

    # For SPlus compatibility:
    if (!exists("beta"))
        beta <- function (a, b) exp( lgamma(a) + lgamma(b) -lgamma(a+b) )

    # Generate Random Deviates:
    weight = xi / (xi + 1/xi)
    z = runif(n, -weight, 1-weight)
    Xi = xi^sign(z)
    Random = -abs(rstd(n, nu = nu))/Xi * sign(z)

    # Scale:
    m1 = 2 * sqrt(nu-2) / (nu-1) / beta(1/2, nu/2)
    mu = m1*(xi-1/xi)
    sigma =  sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)
    Random = (Random - mu ) / sigma

    # Return value:
    Random
}



################################################################################

