
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
#  nigMED               Returns true NIG median
#  nigIQR               Returns true NIG inter quartal range
#  nigSKEW              Returns true NIG robust skewness
#  nigKURT              Returns true NIG robust kurtosis
################################################################################
 

nigMED <-  
function(alpha=1, beta=0, delta=1, mu=0) 
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns true NIG median
    
    # FUNCTION:
    
    # NIG Median:
    Q = .qnigC(p=0.5, alpha, beta, delta, mu)  
    med = c(MED = Q) 
    
    # Return Value:
    med
}
  

# ------------------------------------------------------------------------------
 
    
nigIQR <- 
function(alpha=1, beta=0, delta=1, mu=0) 
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns true NIG inter quartal range
    
    # FUNCTION:
    
    # NIG Inter Quartile Range:
    Q = numeric()
    Q[1] = .qnigC(p=0.25, alpha, beta, delta, mu)  
    Q[2] = .qnigC(p=0.75, alpha, beta, delta, mu) 
    iqr = c(IQR = Q[[2]] - Q[[1]]) 
    
    # Return Value:
    iqr
}
 

# ------------------------------------------------------------------------------
 
  
nigSKEW <- 
function(alpha=1, beta=0, delta=1, mu=0) 
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns true NIG robust skewness
    
    # FUNCTION:
    
    # NIG Robust Skewness:
    Q = numeric()
    Q[1] = .qnigC(p=0.25, alpha, beta, delta, mu)
    Q[2] = .qnigC(p=0.50, alpha, beta, delta, mu) 
    Q[3] = .qnigC(p=0.75, alpha, beta, delta, mu)
    skew = c(SKEW = ( Q[[3]] + Q[[1]] - 2* Q[[2]] ) / (Q[[3]] - Q[[1]] ) ) 
    
    # Return Value:
    skew
}
   

# ------------------------------------------------------------------------------
 

nigKURT <- 
function(alpha=1, beta=0, delta=1, mu=0) 
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns true NIG robust kurtosis
    
    # FUNCTION:
    
    # NIG Robust Kurtosis:
    Q = numeric()
    for (p in (1:7)/8) Q = c(Q, .qnigC(p, alpha, beta, delta, mu))
    kurt = c(KURT = (Q[[7]] - Q[[5]] + Q[[3]] - Q[[1]]) / (Q[[6]] - Q[[2]])) 
    
    # Return Value:
    kurt
}


################################################################################

