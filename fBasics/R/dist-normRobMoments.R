
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
#  normMED               Returns true normal median
#  normIQR               Returns true normal inter quartal range
#  normSKEW              Returns true normal robust skewness
#  normKURT              Returns true normal robust kurtosis
################################################################################
 

normMED <-  
function(mean =0, sd = 1) 
{
    # A function implemented by Diethelm Wuertz 
    
    # Description:
    #   Returns true gh median
    
    # Arguments:
    #   mean - a numeric value, the location parameter
    #   sd - a numeric value, the scale parameter
    
    # FUNCTION:
    
    # Median:
    Q = qnorm(p=0.5, mean, sd)  
    med = c(MED = Q) 
    
    # Return Value:
    med
}
  

# ------------------------------------------------------------------------------
 
    
normIQR <- 
function(mean=0, sd=1) 
{
    # A function implemented by Diethelm Wuertz 
    
    # Description:
    #   Returns true gh inter quartal range
    
    # Arguments:
    #   mean - a numeric value, the location parameter
    #   sd - a numeric value, the scale parameter
 
    # FUNCTION:
    
    # Inter Quartile Range:
    Q = numeric()
    Q[1] = qnorm(p=0.25, mean, sd)  
    Q[2] = qnorm(p=0.75, mean, sd) 
    iqr = c(IQR = Q[[2]] - Q[[1]]) 
    
    # Return Value:
    iqr
}
 

# ------------------------------------------------------------------------------
 
  
normSKEW <- 
function(mean=0, sd=1) 
{
    # A function implemented by Diethelm Wuertz 
    
    # Description:
    #   Returns true gh robust gh skewness
    
    # Arguments:
    #   mean - a numeric value, the location parameter
    #   sd - a numeric value, the scale parameter
    
    # FUNCTION:
    
    # Robust Skewness:
    Q = numeric()
    Q[1] = qnorm(p=0.25, mean, sd)
    Q[2] = qnorm(p=0.50, mean, sd) 
    Q[3] = qnorm(p=0.75, mean, sd)
    skew = c(SKEW = ( Q[[3]] + Q[[1]] - 2* Q[[2]] ) / (Q[[3]] - Q[[1]] ) ) 
    
    # Return Value:
    skew
}
   

# ------------------------------------------------------------------------------
 

normKURT <- 
function(mean=0, sd=1) 
{
    # A function implemented by Diethelm Wuertz 
    
    # Description:
    #   Returns true gh robust gh kurtosis
    
    # Arguments:
    #   mean - a numeric value, the location parameter
    #   sd - a numeric value, the scale parameter
    #   mu - a numeric value, the first shape parameter
    #   nu - a numeric value, the second parameter
    
    # FUNCTION:
    
    # Robust Kurtosis:
    Q = numeric()
    for (p in (1:7)/8) Q = c(Q, qnorm(p, mean, sd))
    kurt = c(KURT = (Q[[7]] - Q[[5]] + Q[[3]] - Q[[1]]) / (Q[[6]] - Q[[2]])) 
    
    # Return Value:
    kurt
}


################################################################################

