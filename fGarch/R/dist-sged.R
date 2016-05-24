
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
# FUNCTION:              DESCRIPTION:
#  dsged                  Density for the skewed GED
#  psged                  Probability function for the skewed GED
#  qsged                  Quantile function for the skewed GED
#  rsged                  Random Number Generator for the skewed GED
# FUNCTION:              DESCRIPTION:
#  .dsged                 Internal, density for the skewed GED
#  .psged                 Internal, probability function for the skewed GED
#  .qsged                 Internal, quantile function for the skewed GED
#  .rsged                 Internal, random Number Generator for the skewed GED
################################################################################

      
dsged <- 
function(x, mean = 0, sd = 1, nu = 2, xi = 1.5, log = FALSE)
{   
    # A function implemented by Diethelm Wuertz 

    # Description:
    #   Compute the density function of the 
    #   skewed generalized error distribution
    
    # FUNCTION:
    
    # Params:
    if (length(mean) == 4) {
        xi = mean[4]
        nu = mean[3]
        sd = mean[2]
        mean = mean[1]
    } 
        
    # Shift and Scale:
    result = .dsged(x = (x-mean)/sd, nu = nu, xi = xi) / sd
    
    # Log:
    if(log) result = log(result)
    
    # Return Value:
    result
}


# ------------------------------------------------------------------------------

      
psged <- 
function(q, mean = 0, sd = 1, nu = 2, xi = 1.5)
{   
    # A function implemented by Diethelm Wuertz 

    # Description:
    #   Compute the distribution function of the 
    #   skewed generalized error distribution
    
    # FUNCTION:
              
    # Shift and Scale:
    result = .psged(q = (q-mean)/sd, nu = nu, xi = xi)
          
    # Return Value:
    result
}


# ------------------------------------------------------------------------------    

        
qsged <- 
function(p, mean = 0, sd = 1, nu = 2, xi = 1.5)
{   
    # A function implemented by Diethelm Wuertz 

    # Description:
    #   Compute the quantile function of the 
    #   skewed generalized error distribution
    
    # FUNCTION:
        
    # Shift and Scale:
    result = .qsged(p = p, nu = nu, xi = xi) * sd + mean
    
    # Return Value:
    result
}

    
# ------------------------------------------------------------------------------


rsged <- 
function(n, mean = 0, sd = 1, nu = 2, xi = 1.5)
{   
    # A function implemented by Diethelm Wuertz 

    # Description:
    #   Generate random deviates from the 
    #   skewed generalized error distribution
    
    # FUNCTION:
        
    # Shift and Scale:
    result = .rsged(n = n, nu = nu, xi = xi) * sd + mean
    
    # Return Value:
    result
}


################################################################################


.dsged <-  
function(x, nu, xi) 
{   
    # A function implemented by Diethelm Wuertz 

    # Description:
    #   Internal Function
    
    # FUNCTION:
    
    # Standardize:
    lambda = sqrt ( 2^(-2/nu) * gamma(1/nu) / gamma(3/nu) )
    g  = nu / ( lambda * (2^(1+1/nu)) * gamma(1/nu) )
    m1 = 2^(1/nu) * lambda * gamma(2/nu) / gamma(1/nu)
    mu = m1*(xi-1/xi)
    sigma =  sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)
    z = x*sigma + mu  
    
    # Compute:
    Xi = xi^sign(z)
    g = 2  / (xi + 1/xi)    
    Density = g * dged(x = z/Xi, nu=nu)  
    
    # Return Value:
    Density * sigma 
}


# ------------------------------------------------------------------------------


.psged <- 
function(q, nu, xi) 
{   
    # A function implemented by Diethelm Wuertz 

    # Description:
    #   Internal Function
    
    # FUNCTION:
    
    # Standardize:
    lambda = sqrt ( 2^(-2/nu) * gamma(1/nu) / gamma(3/nu) )
    g  = nu / ( lambda * (2^(1+1/nu)) * gamma(1/nu) )
    m1 = 2^(1/nu) * lambda * gamma(2/nu) / gamma(1/nu)
    mu = m1*(xi-1/xi)
    sigma =  sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)
    z = q*sigma + mu
    
    # Compute:  
    Xi = xi^sign(z)
    g = 2  / (xi + 1/xi)    
    Probability = Heaviside(z) - sign(z) * g * Xi * pged(q = -abs(z)/Xi, nu=nu)
    
    # Return Value:
    Probability 
}


# ------------------------------------------------------------------------------


.qsged <- 
function(p, nu, xi) 
{   
    # A function implemented by Diethelm Wuertz 

    # Description:
    #   Internal Function
    
    # FUNCTION:
    
    # Standardize:
    lambda = sqrt ( 2^(-2/nu) * gamma(1/nu) / gamma(3/nu) )
    g  = nu / ( lambda * (2^(1+1/nu)) * gamma(1/nu) )
    m1 = 2^(1/nu) * lambda * gamma(2/nu) / gamma(1/nu)
    mu = m1*(xi-1/xi)
    sigma =  sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)
    
    # Compute:  
    g = 2  / (xi + 1/xi)
    sig = sign(p-1/2) 
    Xi = xi^sig       
    p = (Heaviside(p-1/2)-sig*p) / (g*Xi)
    Quantile = (-sig*qged(p=p, sd=Xi, nu=nu) - mu ) / sigma
    
    # Return Value:
    Quantile 
}


# ------------------------------------------------------------------------------


.rsged <- 
function(n, nu, xi) 
{   
    # A function implemented by Diethelm Wuertz 

    # Description:
    #   Internal Function
    
    # FUNCTION:
    
    # Generate Random Deviates:
    weight = xi / (xi + 1/xi)
    z = runif(n, -weight, 1-weight)
    Xi = xi^sign(z)
    Random = -abs(rged(n, nu=nu))/Xi * sign(z)  
    
    # Scale:
    lambda = sqrt ( 2^(-2/nu) * gamma(1/nu) / gamma(3/nu) )
    g  = nu / ( lambda * (2^(1+1/nu)) * gamma(1/nu) )
    m1 = 2^(1/nu) * lambda * gamma(2/nu) / gamma(1/nu)
    mu = m1*(xi-1/xi)
    sigma =  sqrt((1-m1^2)*(xi^2+1/xi^2) + 2*m1^2 - 1)
    Random = (Random - mu ) / sigma 
    
    # Return value:
    Random 
}


################################################################################


