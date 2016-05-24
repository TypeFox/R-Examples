
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
# FUNCTION:            DESCRIPTION:
#  ghtMED               Returns true GHT median
#  ghtIQR               Returns true GHT inter quartal range
#  ghtSKEW              Returns true GHT robust skewness
#  ghtKURT              Returns true GHT robust kurtosis
################################################################################
 

ghtMED <-  
function(beta = 0.1, delta = 1, mu = 0, nu = 10) 
{
    # A function implemented by Diethelm Wuertz 
    
    # Description:
    #   Returns true ght median
    
    # Arguments:
    #   beta - a numeric value, the location parameter
    #   delta - a numeric value, the scale parameter
    #   mu - a numeric value, the first shape parameter
    #   nu - a numeric value, the second parameter
    
    # FUNCTION:
    
    # ght Median
    Q = qght(p=0.5, beta, delta, mu, nu)  
    med = c(MED = Q) 
    
    # Return Value:
    med
}
  

# ------------------------------------------------------------------------------
 
    
ghtIQR <- 
function(beta = 0.1, delta = 1, mu = 0, nu = 10) 
{
    # A function implemented by Diethelm Wuertz 
    
    # Description:
    #   Returns true ght inter quartal range
    
    # Arguments:
    #   beta - a numeric value, the location parameter
    #   delta - a numeric value, the scale parameter
    #   mu - a numeric value, the first shape parameter
    #   nu - a numeric value, the second parameter
    
    # FUNCTION:
    
    # ght Inter Quartile Range
    Q = numeric()
    Q[1] = qght(p=0.25, beta, delta, mu, nu)  
    Q[2] = qght(p=0.75, beta, delta, mu, nu) 
    iqr = c(IQR = Q[[2]] - Q[[1]]) 
    
    # Return Value:
    iqr
}
 

# ------------------------------------------------------------------------------
 
  
ghtSKEW <- 
function(beta = 0.1, delta = 1, mu = 0, nu = 10) 
{
    # A function implemented by Diethelm Wuertz 
    
    # Description:
    #   Returns true ght robust ght skewness
    
    # Arguments:
    #   beta - a numeric value, the location parameter
    #   delta - a numeric value, the scale parameter
    #   mu - a numeric value, the first shape parameter
    #   nu - a numeric value, the second parameter
    
    # FUNCTION:
    
    # ght Robust Skewness:
    Q = numeric()
    Q[1] = qght(p=0.25, beta, delta, mu, nu)
    Q[2] = qght(p=0.50, beta, delta, mu, nu) 
    Q[3] = qght(p=0.75, beta, delta, mu, nu)
    skew = c(SKEW = ( Q[[3]] + Q[[1]] - 2* Q[[2]] ) / (Q[[3]] - Q[[1]] ) ) 
    
    # Return Value:
    skew
}
   

# ------------------------------------------------------------------------------
 

ghtKURT <- 
function(beta = 0.1, delta = 1, mu = 0, nu = 10) 
{
    # A function implemented by Diethelm Wuertz 
    
    # Description:
    #   Returns true ght robust ght kurtosis
    
    # Arguments:
    #   beta - a numeric value, the location parameter
    #   delta - a numeric value, the scale parameter
    #   mu - a numeric value, the first shape parameter
    #   nu - a numeric value, the second parameter
    
    # FUNCTION:
    
    # ght Robust Kurtosis:
    Q = numeric()
    for (p in (1:7)/8) Q = c(Q, qght(p, beta, delta, mu, nu))
    kurt = c(KURT = (Q[[7]] - Q[[5]] + Q[[3]] - Q[[1]]) / (Q[[6]] - Q[[2]])) 
    
    # Return Value:
    kurt
}


################################################################################

