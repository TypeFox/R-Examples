
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
# FUNCTION:             BIVARIATE CAUCHY DISTRIBUTION:
#  pcauchy2d             Computes bivariate Cauchy probability function
#  dcauchy2d             Computes bivariate Cauchy density function
#  rcauchy2d             Generates bivariate Cauchy random deviates
################################################################################


pcauchy2d = 
function(x, y = x, rho = 0) 
{   # A function Implemented by Diethelm Wuertz

    # Description:
    #   Computes bivariate Cauchy probability function
    
    # Arguments:
    #   x, y - two numeric values or vectors of the same length at
    #       which the probability will be computed. 
    
    # Example:
    #   pt2d(rnorm(5), rnorm(5), 0.5, 5)
    
    # Value:
    #   returns a numeric vector of probabilities of the same length
    #   as the input vectors
   
    # FUNCTION:
    
    # Settings:
    # Probaility:
    ans  = pt2d(x = x, y = y, rho = rho, nu = 1) 
    attr(ans, "control") = c(rho = rho)
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


dcauchy2d = 
function(x, y = x, rho = 0)
{   # A function implemented by Diethelm Wuertz

    # Arguments:
    #   n - number of random deviates to be generated
    #   rho - the linear correlation, a numeric value between 
    #       minus one and one.
    
    # Description:
    #   Computes bivariate Cauchy density function
    
    # Note:
    #   Partly copied from contributed R package 'sn'
    
    # FUNCTION:
    
    # Density:
    density = dt2d(x = x, y = y, rho = rho, nu = 1)
    attr(density, "control") = c(rho = rho)
    
    # Return value:
    density
}


# ------------------------------------------------------------------------------


rcauchy2d =
function(n, rho = 0) 
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Generates bivariate Cauchy random deviates
    
    # Arguments:
    #   n - number of random deviates to be generated
    #   rho - the linear correlation, a numeric value between 
    #       minus one and one.
    
    # Note:
    #   Partly copied from contributed R package 'mvtnorm'
    #   Author Friedrich Leisch
    
    # FUNCTION:
    
    # Random Deviates:
    ans = rt2d(n = n, rho = rho)
    attr(ans, "control") = c(rho = rho)
    
    # Return Value:
    ans
}


################################################################################

