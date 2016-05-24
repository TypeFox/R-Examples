
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
#  gldMode               Computes the Generalized Lambda Distribution mode
################################################################################


gldMode <- 
function(lambda1 = 0, lambda2 = -1, lambda3 = -1/8, lambda4 = -1/8)
{   
    # A function implemented by Diethelm Wuertz
    
    # Arguments:
    #   lambda1 - a numeric value, the location parameter
    #   lambda2 - a numeric value, the scale parameter
    #   lambda3 - a numeric value, the first shape parameter
    #   lambda4 - a numeric value, the second parameter
    
    # Description:
    #   Computes the Generalized Lambda Distribution mode
    
    # FUNCTION:
    
    # Find Maximum: 
    min = qgld(0.01, lambda1, lambda2, lambda3, lambda4)
    max = qgld(0.99, lambda1, lambda2, lambda3, lambda4)
    ans = optimize(f = dgld, interval = c(min, max), 
        lambda1 = lambda1, lambda2 = lambda2, 
        lambda3 = lambda3, lambda4 = lambda4,  
        maximum = TRUE, tol = .Machine$double.eps)$maximum
    
    # Return Value:
    ans
}


################################################################################

