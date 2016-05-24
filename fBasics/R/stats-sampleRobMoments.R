
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
# FUNCTION:             DESCRIPTION:
#  sampleMED             Returns sample median
#  sampleIQR             Returns sample inter quartal range
#  sampleSKEW            Returns robust sample skewness
#  sampleKURT            Returns robust sample kurtosis
################################################################################


sampleMED <-
function(x)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns sample median
    
    # FUNCTION:
    
    # Sample Median:
    x = as.vector(x)
    med = c(MED = quantile(x, 0.50)[[1]])
    
    # Return Value:
    med
}


# ------------------------------------------------------------------------------


sampleIQR = 
function(x) 
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns sample inter quartal range
    
    # FUNCTION:
    
    # Sample Inter Quartile Range
    x = as.vector(x)
    iqr = c(IQR = diff(quantile(x, c(0.25, 0.75)))[[1]])
    
    # Return Value:
    iqr
}


# ------------------------------------------------------------------------------


sampleSKEW = 
function(x) 
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns robust sample skewness
    
    # FUNCTION:
    
    # Robust Sample Skew
    x = as.vector(x)
    q = quantile(x, c(0.25, 0.5, 0.75))
    skew = c(SKEW = q[[3]] + q[[1]] - 2*q[[2]]) / (q[[3]]-q[[1]]) 
    
    # Return Value:
    skew
}


# ------------------------------------------------------------------------------


sampleKURT = 
function(x) 
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns robust sample kurtosis
    # FUNCTION:
    
    # Robust Sample Kurtosis
    x = as.vector(x)
    q = quantile(x, (1:7)/8)
    kurt = c(KURT = (q[[7]] - q[[5]] + q[[3]] - q[[1]]) / (q[[6]] - q[[2]])) 
    
    # Return Value:
    kurt
}


################################################################################
    
  