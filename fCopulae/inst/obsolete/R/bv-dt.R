
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
# FUNCTION:             BIVARIATE STUDENT-T DISTRIBUTION:
#  pt2d                  Computes bivariate Student-t probability function
#  dt2d                  Computes bivariate Student-t density function
#  rt2d                  Generates bivariate Student-t random deviates
################################################################################


pt2d = 
function(x, y = x, rho = 0, nu = 4) 
{   # pnorm2d: A copy from R package "sn"

    # Description:
    #   Computes bivariate Student-t probability function
    
    # Arguments:
    #   x, y - two numeric values or vectors of the same length at
    #       which the probability will be computed. 
    
    # Example:
    #   pt2d(rnorm(5), rnorm(5), 0.5, 5)
    
    # Value:
    #   returns a numeric vector of probabilities of the same length
    #   as the input vectors
   
    # FUNCTION:
    
    # Normal Limit:
    if (nu == Inf) return(pnorm2d(x = x, y = y, rho = rho)) 
    
    # Settings:
    sigma = diag(2) 
    sigma[1, 2] = sigma[2, 1] = rho 
    X = cbind(x, y)
    
    # Probaility:
    ans  = pmvst(X, dim = 2, mu = c(0, 0), Omega = sigma, 
        alpha = c(0, 0), df = nu) 
    attr(ans, "control") = c(rho = rho, nu = nu)
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


dt2d = 
function(x, y = x, rho = 0, nu = 4)
{   # A function implemented by Diethelm Wuertz

    # Arguments:
    #   n - number of random deviates to be generated
    #   rho - the linear correlation, a numeric value between 
    #       minus one and one.
    
    # Description:
    #   Computes bivariate Student-t density function
    
    # Example:
    #   dt2d(rnorm(5), rnorm(5), 0.5, 5)
    
    # Note:
    #   Partly copied from contributed R package 'sn'
    
    # FUNCTION:
    
    # Normal Limit:
    if (nu == Inf) return(dnorm2d(x = x, y = y, rho = rho)) 
    
    # Argument:
    xoy = (x^2 - 2*rho*x*y + y^2)/ (2*(1 - rho^2))
    
    # Density:
    density = (1 + 2*xoy/nu)^(-(nu+2)/2) / (2*pi*sqrt(1-rho^2))
    attr(density, "control") = c(rho = rho, nu = nu)
    
    # Return value:
    density
}


# ------------------------------------------------------------------------------


rt2d =
function(n, rho = 0, nu = 4) 
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Generates bivariate Student-t random deviates
    
    # Arguments:
    #   n - number of random deviates to be generated
    #   rho - the linear correlation, a numeric value between 
    #       minus one and one.
    
    # Note:
    #   Partly copied from contributed R package 'mvtnorm'
    #   Author Friedrich Leisch
    
    # FUNCTION:
    
    # Normal Limit:
    if (nu == Inf) return(rnorm2d(n = n, rho = rho)) 
    
    # Random Deviates:
    ans = rnorm2d(n, rho)/sqrt(rchisq(n, nu)/nu)
    attr(ans, "control") = c(rho = rho, nu = nu)
    
    # Return Value:
    ans
}


################################################################################

