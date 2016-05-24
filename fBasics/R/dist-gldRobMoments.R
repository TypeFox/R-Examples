
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
#  gldMED               Returns true GLD median
#  gldIQR               Returns true GLD inter quartal range
#  gldSKEW              Returns true GLD robust skewness
#  gldKURT              Returns true GLD robust kurtosis
################################################################################
 

gldMED <-  
function(lambda1=0, lambda2=-1, lambda3=-1/8, lambda4=-1/8) 
{
    # A function implemented by Diethelm Wuertz 
    
    # Description:
    #   Returns true gld median
    
    # Arguments:
    #   lambda1 - a numeric value, the location parameter
    #   lambda2 - a numeric value, the scale parameter
    #   lambda3 - a numeric value, the first shape parameter
    #   lambda4 - a numeric value, the second parameter
    
    # FUNCTION:
    
    # gld Median
    Q = qgld(p=0.5, lambda1, lambda2, lambda3, lambda4)  
    med = c(MED = Q) 
    
    # Return Value:
    med
}
  

# ------------------------------------------------------------------------------
 
    
gldIQR <- 
function(lambda1=0, lambda2=-1, lambda3=-1/8, lambda4=-1/8) 
{
    # A function implemented by Diethelm Wuertz 
    
    # Description:
    #   Returns true gld inter quartal range
    
    # Arguments:
    #   lambda1 - a numeric value, the location parameter
    #   lambda2 - a numeric value, the scale parameter
    #   lambda3 - a numeric value, the first shape parameter
    #   lambda4 - a numeric value, the second parameter
    
    # FUNCTION:
    
    # gld Inter Quartile Range
    Q = numeric()
    Q[1] = qgld(p=0.25, lambda1, lambda2, lambda3, lambda4)  
    Q[2] = qgld(p=0.75, lambda1, lambda2, lambda3, lambda4) 
    iqr = c(IQR = Q[[2]] - Q[[1]]) 
    
    # Return Value:
    iqr
}
 

# ------------------------------------------------------------------------------
 
  
gldSKEW <- 
function(lambda1=0, lambda2=-1, lambda3=-1/8, lambda4=-1/8) 
{
    # A function implemented by Diethelm Wuertz 
    
    # Description:
    #   Returns true gld robust gld skewness
    
    # Arguments:
    #   lambda1 - a numeric value, the location parameter
    #   lambda2 - a numeric value, the scale parameter
    #   lambda3 - a numeric value, the first shape parameter
    #   lambda4 - a numeric value, the second parameter
    
    # FUNCTION:
    
    # gld Robust Skewness:
    Q = numeric()
    Q[1] = qgld(p=0.25, lambda1, lambda2, lambda3, lambda4)
    Q[2] = qgld(p=0.50, lambda1, lambda2, lambda3, lambda4) 
    Q[3] = qgld(p=0.75, lambda1, lambda2, lambda3, lambda4)
    skew = c(SKEW = ( Q[[3]] + Q[[1]] - 2* Q[[2]] ) / (Q[[3]] - Q[[1]] ) ) 
    
    # Return Value:
    skew
}
   

# ------------------------------------------------------------------------------
 

gldKURT <- 
function(lambda1=0, lambda2=-1, lambda3=-1/8, lambda4=-1/8) 
{
    # A function implemented by Diethelm Wuertz 
    
    # Description:
    #   Returns true gld robust gld kurtosis
    
    # Arguments:
    #   lambda1 - a numeric value, the location parameter
    #   lambda2 - a numeric value, the scale parameter
    #   lambda3 - a numeric value, the first shape parameter
    #   lambda4 - a numeric value, the second parameter
    
    # FUNCTION:
    
    # gld Robust Kurtosis:
    Q = numeric()
    for (p in (1:7)/8) Q = c(Q, qgld(p, lambda1, lambda2, lambda3, lambda4))
    kurt = c(KURT = (Q[[7]] - Q[[5]] + Q[[3]] - Q[[1]]) / (Q[[6]] - Q[[2]])) 
    
    # Return Value:
    kurt
}


################################################################################

