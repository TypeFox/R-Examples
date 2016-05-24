
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
# FUNCTION:           DESCRIPTION:
#  ghMED               Returns true gh median
#  ghIQR               Returns true gh inter quartal range
#  ghSKEW              Returns true gh robust skewness
#  ghKURT              Returns true gh robust kurtosis
################################################################################
 

ghMED <-  
function(alpha = 1, beta = 0, delta = 1, mu = 0, lambda = -1/2) 
{
    # A function implemented by Diethelm Wuertz 
    
    # Description:
    #   Returns true gh median
    
    # Arguments:
    #   beta - a numeric value, the location parameter
    #   delta - a numeric value, the scale parameter
    #   mu - a numeric value, the first shape parameter
    #   nu - a numeric value, the second parameter
    
    # FUNCTION:
    
    # gh Median
    Q = qgh(p=0.5, alpha, beta, delta, mu, lambda)  
    med = c(MED = Q) 
    
    # Return Value:
    med
}
  

# ------------------------------------------------------------------------------
 
    
ghIQR <- 
function(alpha = 1, beta = 0, delta = 1, mu = 0, lambda = -1/2) 
{
    # A function implemented by Diethelm Wuertz 
    
    # Description:
    #   Returns true gh inter quartal range
    
    # Arguments:
    #   beta - a numeric value, the location parameter
    #   delta - a numeric value, the scale parameter
    #   mu - a numeric value, the first shape parameter
    #   nu - a numeric value, the second parameter
    
    # FUNCTION:
    
    # gh Inter Quartile Range
    Q = numeric()
    Q[1] = qgh(p=0.25, alpha, beta, delta, mu, lambda)  
    Q[2] = qgh(p=0.75, alpha, beta, delta, mu, lambda) 
    iqr = c(IQR = Q[[2]] - Q[[1]]) 
    
    # Return Value:
    iqr
}
 

# ------------------------------------------------------------------------------
 
  
ghSKEW <- 
function(alpha = 1, beta = 0, delta = 1, mu = 0, lambda = -1/2) 
{
    # A function implemented by Diethelm Wuertz 
    
    # Description:
    #   Returns true gh robust gh skewness
    
    # Arguments:
    #   beta - a numeric value, the location parameter
    #   delta - a numeric value, the scale parameter
    #   mu - a numeric value, the first shape parameter
    #   nu - a numeric value, the second parameter
    
    # FUNCTION:
    
    # gh Robust Skewness:
    Q = numeric()
    Q[1] = qgh(p=0.25, alpha, beta, delta, mu, lambda)
    Q[2] = qgh(p=0.50, alpha, beta, delta, mu, lambda) 
    Q[3] = qgh(p=0.75, alpha, beta, delta, mu, lambda)
    skew = c(SKEW = ( Q[[3]] + Q[[1]] - 2* Q[[2]] ) / (Q[[3]] - Q[[1]] ) ) 
    
    # Return Value:
    skew
}
   

# ------------------------------------------------------------------------------
 

ghKURT <- 
function(alpha = 1, beta = 0, delta = 1, mu = 0, lambda = -1/2) 
{
    # A function implemented by Diethelm Wuertz 
    
    # Description:
    #   Returns true gh robust gh kurtosis
    
    # Arguments:
    #   beta - a numeric value, the location parameter
    #   delta - a numeric value, the scale parameter
    #   mu - a numeric value, the first shape parameter
    #   nu - a numeric value, the second parameter
    
    # FUNCTION:
    
    # gh Robust Kurtosis:
    Q = numeric()
    for (p in (1:7)/8) Q = c(Q, qgh(p, alpha, beta, delta, mu, lambda))
    kurt = c(KURT = (Q[[7]] - Q[[5]] + Q[[3]] - Q[[1]]) / (Q[[6]] - Q[[2]])) 
    
    # Return Value:
    kurt
}


################################################################################

