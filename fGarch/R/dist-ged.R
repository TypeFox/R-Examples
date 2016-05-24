
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
# FUNCTION:              GED DISTRIBUTION:
#  dged                   Density for the Generalized Error Distribution
#  pged                   Probability function for the GED
#  qged                   Quantile function for the GED
#  rged                   Random Number Generator for the GED
################################################################################


dged <- 
function(x, mean = 0, sd = 1, nu = 2, log = FALSE)
{   
    # A function imlemented by Diethelm Wuertz

    # Description:
    #   Compute the density for the 
    #   generalized error distribution.
    
    # FUNCTION:
    
    # Params:
    if (length(mean) == 3) {
        nu = mean[3]
        sd = mean[2]
        mean = mean[1]
    }   
    
    # Compute Density:
    z = (x - mean ) / sd
    lambda = sqrt ( 2^(-2/nu) * gamma(1/nu) / gamma(3/nu) )
    g  = nu / ( lambda * (2^(1+1/nu)) * gamma(1/nu) )
    result = g * exp (-0.5*(abs(z/lambda))^nu) / sd
    
    # Log:
    if(log) result = log(result)
    
    # Return Value
    result
}

        
# ------------------------------------------------------------------------------


pged <-  
function(q, mean = 0, sd = 1, nu = 2)
{   
    # A function implemented by Diethelm Wuertz
        
    # Description:
    #   Compute the probability for the  
    #   generalized error distribution.
    
    # FUNCTION:
        
    # Compute Probability:
    q = (q - mean ) / sd
    lambda = sqrt ( 2^(-2/nu) * gamma(1/nu) / gamma(3/nu) )
    g  = nu / ( lambda * (2^(1+1/nu)) * gamma(1/nu) )
    h = 2^(1/nu) * lambda * g * gamma(1/nu) / nu
    s = 0.5 * ( abs(q) / lambda )^nu
    result = 0.5 + sign(q) * h * pgamma(s, 1/nu)
    
    # Return Value:
    result
}
        

# ------------------------------------------------------------------------------


qged <- 
function(p, mean = 0, sd = 1, nu = 2)
{   
    # A function implemented by Diethelm Wuertz
        
    # Description:
    #   Compute the quantiles for the  
    #   generalized error distribution.
    
    # FUNCTION:
    
    # Compute Quantiles:
    lambda = sqrt ( 2^(-2/nu) * gamma(1/nu) / gamma(3/nu) )
    q = lambda * (2*qgamma((abs(2*p-1)), 1/nu))^(1/nu)
    result = q*sign(2*p-1) * sd + mean
    
    # Return Value:
    result
}


# ------------------------------------------------------------------------------

    
rged <-  
function(n, mean = 0, sd = 1, nu = 2)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Generate GED random deviates. The function uses the 
    #   method based on the transformation of a Gamma random 
    #   variable.
    
    # FUNCTION:
    
    # Generate Random Deviates:
    lambda = sqrt ( 2^(-2/nu) * gamma(1/nu) / gamma(3/nu) )
    # print(lambda)
    r = rgamma(n, 1/nu)
    z =  lambda * (2*r)^(1/nu) * sign(runif(n)-1/2)
    result = z * sd + mean

    
    # Return Value:
    result
}


################################################################################

