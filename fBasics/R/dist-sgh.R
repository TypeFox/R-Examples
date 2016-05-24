
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
# You should have received A copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:            DESCRIPTION:
#  dsgh                 Returns density of the standardized gh distribution
#  psgh                 Returns probabilities of the sgh distribution
#  qsgh                 Returns quantiles of the sgh distribution
#  rsgh                 Generates sgh distributed random variates
# FUNCTION:            DESCRIPTION:
#  .kappaGH             Returns modified Bessel function ratio
#  .deltaKappaGH        Returns difference of Bessel function ratios
#  .paramGH             Change parameterization to alpha, beta, delta mu
################################################################################


dsgh <-  
function(x, zeta = 1, rho = 0, lambda = 1, log = FALSE) 
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns density of the sgh distribution
    
    # FUNCTION:    
    
    # Parameters:
    if (length(zeta) == 3) {
       lambda = zeta[3]
       rho = zeta[2]
       zeta = zeta[1]
    } 
    
    # Compute Density:
    param = .paramGH(zeta, rho, lambda)
    ans = dgh(x, param[1], param[2], param[3], param[4], lambda, log)
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


psgh <-  
function(q, zeta = 1, rho = 0, lambda = 1) 
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns probabilities of the sgh distribution
    
    # FUNCTION:
    
    # Compute Probabilities:
    param = .paramGH(zeta, rho, lambda)
    ans = pgh(q, param[1], param[2], param[3], param[4], lambda)
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


qsgh <-  
function(p, zeta = 1, rho = 0, lambda = 1) 
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns quantiles of the sgh distribution
    
    # FUNCTION:
      
    # Compute Quantiles:
    param = .paramGH(zeta, rho, lambda)
    ans = qgh(p, param[1], param[2], param[3], param[4], lambda)
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


rsgh <-  
function(n, zeta = 1, rho = 0, lambda = 1) 
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Generates sgh distributed random variates
    
    # FUNCTION:
    
    # Generate Random Numbers:
    param = .paramGH(zeta, rho, lambda)
    ans = rgh(n, param[1], param[2], param[3], param[4], lambda)
    
    # Return Value:
    ans
}


################################################################################


.kappaGH <- 
function(x, lambda = 1)
{    
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns modified Bessel function ratio
    
    # FUNCTION:
    
    # Check:
    stopifnot(x >= 0)
    stopifnot(length(lambda) == 1)
    
    # Ratio:
    if (lambda == -0.5) {
        # NIG:
        kappa = 1/x
    } else {
        # GH:
        kappa = (
            besselK(x, lambda+1, expon.scaled = TRUE) /
            besselK(x, lambda, expon.scaled = TRUE) ) / x
    }
    
    # Return Value:
    kappa
}


# ------------------------------------------------------------------------------


.deltaKappaGH <-  
function(x, lambda = 1)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns difference of Bessel functions ratios
    
    # FUNCTION:
    
    # Difference in Ratios:
    if (lambda == -0.5) {
        # NIG:
        # Replace this with the recursion relation ...
        deltaKappa = .kappaGH(x, lambda+1) - .kappaGH(x, lambda)
    } else {
        # GH:
        deltaKappa = .kappaGH(x, lambda+1) - .kappaGH(x, lambda)
    }
    
    # Return Value:
    deltaKappa
}


################################################################################


.paramGH <-  
function(zeta = 1, rho = 0 , lambda = 1)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Change parameterizations to alpha(zeta, rho, lambda)
    
    # FUNCTION:
    
    # Transformation:
    Rho2 = 1 - rho^2
    alpha = zeta^2 * .kappaGH(zeta, lambda) / Rho2 
    alpha = alpha * ( 1 + rho^2 * zeta^2 * .deltaKappaGH(zeta, lambda) / Rho2)
    alpha = sqrt(alpha)  
    beta = alpha * rho
    delta = zeta / ( alpha * sqrt(Rho2) )
    mu = -beta * delta^2 * .kappaGH(zeta, lambda)
               
    # Return Value:
    c(alpha = alpha, beta = beta, delta = delta, mu = mu)  
}


################################################################################

