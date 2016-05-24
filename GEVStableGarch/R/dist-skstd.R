
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


##################################################################################
#  FUNCTION:             Skewed version (Fernandez, C. and Steel, M. F. J. (1998))
#                        of the standard t-Student (std) distribution
#
#  dskstd                Density for the skstd
#  pskstd                Probability function for the skstd
#  qskstd                Quantile function for the skstd
#  rskstd                Random Number Generator for the skstd
##################################################################################


dskstd <- 
  function(x, mean = 0, sd = 1, nu = 3, xi = 1, log = FALSE)
  {   
    # A function implemented by Thiago Sousa
    
    # DESCRIPTION:
    #   Compute the density for the 
    #   skewed version of the sandard t-Student.
    #   The skewed t-Student distributio is already defined on
    #   the skewt R package, but it uses th non-standardized version
    #   of the t-Student distribution, i.e., the variance is df/(df-2)
    #   and it does not allow the user to specify the location or 
    #   scale parameter. We are defining this slightly diffeernt 
    #   version just because it is widelly used in the GARCH context.
    #   Reference: Fernandez, C. and Steel, M. F. J. (1998). 
    #   On Bayesian modeling of fat tails and skewness, 
    # J. Am. Statist. Assoc. 93, 359-371.
    
    # INPUT PARAMETERS:
    #  x - the points in which the function will be evaluated
    #  mean - the location parameter
    #  sd - the scale parameter 
    #  nu - the shape parameter ( degrees of freedom) ( > 0)
    #  xi - the asymmetry parameter ( > 0 ). xi = 1 gives a symmetric distr.
    #  log - logical indicating whether or not to report the log of the density
    
    # FUNCTION:
    
    # Params:
    if (length(mean) == 4) {
      xi = mean[4]
      nu = mean[3]
      sd = mean[2]
      mean = mean[1]
    }    
    
    # Error treatment of input parameters
    if(sd <= 0  || nu <= 2 || xi <= 0)
      stop("Failed to verify condition:
           sd <= 0  || nu <= 0 || xi <= 0")
    
    # Compute density points:
    z = (x - mean ) / sd
    result = sqrt( nu / ( nu - 2 ) ) *
             dskt(z * sqrt( nu / ( nu - 2 ) ), df = nu, gamma = xi)
    
    # Log:
    if(log) result = log(result)
    
    # Return Value
    result
  }



# ------------------------------------------------------------------------------


pskstd <- 
  function(q, mean = 0, sd = 1, nu = 3, xi = 1)
  {   
    # A function implemented by Thiago Sousa
    
    # DESCRIPTION:
    #   Compute the probability function for the 
    #   skewed version of the sandard t-Student.
    
    # INPUT PARAMETERS:
    #  q - the points in which the function will be evaluated
    #  mean - the location parameter
    #  sd - the scale parameter 
    #  nu - the shape parameter ( degrees of freedom) ( > 0)
    #  xi - the asymmetry parameter ( > 0 ). xi = 1 gives a symmetric distr.
    #  log - logical indicating whether or not to report the log of the density
    
    # FUNCTION:
    
    # Params:
    if (length(mean) == 4) {
      xi = mean[4]
      nu = mean[3]
      sd = mean[2]
      mean = mean[1]
    }    
    
    # Error treatment of input parameters
    if(sd <= 0  || nu <= 2 || xi <= 0)
      stop("Failed to verify condition:
           sd <= 0  || nu <= 0 || xi <= 0")
    
    # Compute probability points:
    z = (q - mean ) / sd
    result = pskt( z * sqrt( nu / ( nu - 2 ) ), df = nu, gamma = xi)
    
    # Return Value
    result
  }


# ------------------------------------------------------------------------------
  
  
qskstd <- 
  function(p, mean = 0, sd = 1, nu = 3, xi = 1)
  {   
    # A function implemented by Thiago Sousa
    
    # DESCRIPTION:
    #   Compute the quantile function for the 
    #   skewed version of the sandard t-Student.
    
    # INPUT PARAMETERS:
    #  p - the probabilities in which the function will be evaluated
    #  mean - the location parameter
    #  sd - the scale parameter 
    #  nu - the shape parameter ( degrees of freedom) ( > 0)
    #  xi - the asymmetry parameter ( > 0 ). xi = 1 gives a symmetric distr.
    #  log - logical indicating whether or not to report the log of the density
    
    # FUNCTION:
    
    # Params:
    if (length(mean) == 4) {
      xi = mean[4]
      nu = mean[3]
      sd = mean[2]
      mean = mean[1]
    }    
    
    # Error treatment of input parameters
    if(sd <= 0  || nu <= 2 || xi <= 0)
      stop("Failed to verify condition:
           sd <= 0  || nu <= 0 || xi <= 0")
    
    # Compute probability points:
    result = qskt( p, df = nu, gamma = xi) /  sqrt( nu / ( nu - 2 ) )
    
    # Return Value:
    result * sd + mean
  }  

  
  
# ------------------------------------------------------------------------------


rskstd <- 
  function(n, mean = 0, sd = 1, nu = 3, xi = 1)
  {   
    
    # Description:
    #   Generate skstd Random values
    #   using the inverse of the distribution function.
    
    # FUNCTION:
    
    randomUnif = runif(n = n, min = 0, max = 1)
    result = qskstd(p = randomUnif, mean = mean, sd = sd, nu = nu, xi = xi)
    
    # Return Value:
    result
  }


################################################################################


