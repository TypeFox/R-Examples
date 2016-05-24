
# Copyrights (C) 2014 Thiago do Rego Sousa <thiagoestatistico@gmail.com>

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


################################################################################
# FUNCTION:              t3 (now called GAt) distribution proposed in Paolella (1997)
#  pgat                   Probability function for the GAt
#  dgat                   Density for the GAt-distribution 
#  qgat                   Quantile function for the GAt
#  rgat                   Random Number Generator for the GAt
################################################################################


dgat <- 
  function(x, mean = 0, sd = 1, nu = 2, d = 3, xi = 1, log = FALSE)
  {   
    # A function implemented by Thiago Sousa
    
    # Description:
    #   Compute the density for the 
    #   so called t3-distribution (now called GAt) defined in Paolella (1997).
    #   Reference: Paolella M 0886. Tail Estimation and Conditional Modeling 
    #   of Heteroscedastic Time!Series. PhD thesis.
    #   Institute of Statistics and Econometrics. Christian Albrechts University at Kiel
    #   Parameters: mean in R; sd > 0; nu > 0; d > 0; xi > 0; 
    
    # FUNCTION:
    
    # Params:
    if (length(mean) == 5) {
      xi = mean[5]
      d  = mean[4]
      nu = mean[3]
      sd = mean[2]
      mean = mean[1]
    }    
    
    # Error treatment of input parameters
    if(sd <= 0  || nu <= 0 || xi <= 0 || d <= 0)
      stop("Failed to verify condition:
           sd <= 0 || nu <= 0 || xi <= 0 || d <= 0")
    
    # Compute auxiliary variables:
    z = (x - mean ) / sd
    n = length(z)
    arg = z
    indexLessThanZero = which (z < 0, arr.ind = TRUE)
    sizeIndex = length(indexLessThanZero)
    
    # Compute the density points according to their sign
    # all b coefficients are >= 0
    if(sizeIndex == 0) {
      arg = arg / xi  
    # all b coefficients are < 0
    } else if (sizeIndex == n) {
        arg = -arg * xi 
    # default case. we have both pos. and neg. values
    } else if (TRUE) { 
        arg[indexLessThanZero] = -arg[indexLessThanZero] * xi
        arg[-indexLessThanZero] = arg[-indexLessThanZero] / xi 
    }
    
    # Compute density points
    k = ( ( xi + 1/xi ) * 1/d * nu^(1/d) * beta (1/d,nu) )^(-1)
    result = ( k * (1 + (arg^d) / nu )^( -nu-1/d) ) / sd
    # Log:
    if(log) result = log(result)
    
    # Return Value
    result
  }





# ------------------------------------------------------------------------------

pgat <- 
  function(q, mean = 0, sd = 1, nu = 2, d = 3, xi = 1)
  {   
    # A function imlemented by Thiago Sousa
    
    # Description:
    #   Compute the distribution for the 
    #   so called t3-distribution (now called GAt)
    #   Parameters: mean in R; sd > 0; nu > 0; d > 0; xi > 0; 
    
    # FUNCTION:
    
    # Params:
    if (length(mean) == 5) {
      xi = mean[5]
      d  = mean[4]
      nu = mean[3]
      sd = mean[2]
      mean = mean[1]
    }    
    
    # Error treatment of input parameters
    if(sd <= 0  || nu <= 0 || xi <= 0 || d <= 0)
      stop("Failed to verify condition:
           sd <= 0 || nu <= 0 || xi <= 0 || d <= 0")
    
    # Define auxiliary functions
    L <- function (z, nu = nu, d = d, xi = xi)
    {
        nu / ( nu + (-z*xi)^d ) # z must be negative ( <= 0 ), but we do not check it here
    }
    U <- function (z, nu = nu, d = d, xi = xi)
    {
        pw = (z/xi)^d  # z must be negative ( <= 0 ), but we do not check it here
        return ( replace(pw / ( nu + pw ),which (pw == Inf, arr.ind = TRUE),1) )
    } 

    # Compute auxiliary variables:
    z = (q - mean ) / sd
    n = length(z)
    arg = z
    indexLessThanZero = which (z <= 0, arr.ind = TRUE)
    sizeIndex = length(indexLessThanZero)
    
    # Compute distribution points according to their sign
    if(sizeIndex == 0) {
        arg = 1/(1 + xi^2 ) + 1/(1 + xi^(-2) ) * 
            pbeta ( U (z = arg, nu = nu, d = d, xi = xi), 1/d, nu)  
    } else if (sizeIndex == n) {
        arg = 1/(1 + xi^2 ) * 
            pbeta ( L (z = arg, nu = nu, d = d, xi = xi), nu, 1/d)
    } else if (TRUE) { 
        arg[indexLessThanZero] = 1/(1 + xi^2 ) * 
            pbeta ( L (z = arg[indexLessThanZero], nu = nu, d = d, xi = xi), nu, 1/d)
        arg[-indexLessThanZero] = 1/(1 + xi^2 ) + 1/(1 + xi^(-2) ) * 
            pbeta ( U (z = arg[-indexLessThanZero], nu = nu, d = d, xi = xi), 1/d, nu) 
    }
    
    # Return Value
    arg
  }


#------------------------------------------------------------------------------

  
qgat <- 
  function(p, mean = 0, sd = 1, nu = 2, d = 3, xi = 1)  
  {   
    
    # Description:
    #   Compute the quantiles for the  
    #   generalized error distribution using the 
    #   formula for the distribution function and
    #   the quantile function of the Beta Distribution
    #   already available in R
    
    # FUNCTION:
    
    # Define auxiliary functions
    Lp <- function (p = p, nu = nu, d = d, xi = xi)
    {
        qbeta( ( 1 + xi^2 ) * p, nu, 1/d)
    }
    Up <- function (p = p, nu = nu, d = d, xi = xi)
    {
        qbeta( ( p - 1 / ( 1 + xi^2 ) ) * ( 1 + xi^(-2) ), 1/d, nu)
    }
    
    # Compute quantiles located at (-Inf,0] and at (0,+Inf)
    F0 = pgat(0, mean = 0, sd = 1, nu = nu, d = d, xi = xi)
    n = length(p)
    result = rep(NA,n)
    indexLessThanF0 = which (p <= F0, arr.ind = TRUE)
    sizeIndex = length(indexLessThanF0)
    if(sizeIndex == 0) {
      
        U = Up (p = p, nu = nu, d = d, xi = xi)
        result = ( U * nu / (1 - U) )^( 1/d ) * xi
        
    } else if (sizeIndex == n) {
      
         L = Lp (p = p, nu = nu, d = d, xi = xi)
        result = - ( nu/L - nu )^( 1/d ) * 1/xi
        
    } else if (TRUE) {
      
        L = Lp (p = p[indexLessThanF0], nu = nu, d = d, xi = xi)
        U = Up (p = p[-indexLessThanF0], nu = nu, d = d, xi = xi)
        result[indexLessThanF0] = - ( nu/L - nu )^( 1/d ) * 1/xi
        result[-indexLessThanF0] = ( U * nu / (1 - U) )^( 1/d ) * xi 
    }
    
    # Return Value:
    result * sd + mean
  }







# ------------------------------------------------------------------------------


rgat <-  
  function(n, mean = 0, sd = 1, nu = 2, d = 3, xi = 1)  
  {   
    
    # Description:
    #   Generate GAt Random values
    #   using the inverse of the distribution function.
    
    # FUNCTION:
    
    randomUnif = runif(n = n, min = 0, max = 1)
    result = qgat(p = randomUnif, mean = mean, sd = sd, nu = nu, d = d, xi = xi)
    
    # Return Value:
    result
  }


################################################################################

