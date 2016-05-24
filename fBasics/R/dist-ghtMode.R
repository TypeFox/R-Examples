
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
#  ghtMode               Computes the generalized hyperbolic Student-t mode
################################################################################


ghtMode <- 
function(beta = 0.1, delta = 1, mu = 0, nu = 10)
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Computes the mode of the Generalized Hyperbolic Student-t PDF
  
    # Arguments:
    #   beta - skewness parameter
    #   delta - scale parameter
    #   mu - location parameter
    #   nu - shape parameter
     
    # Example:
    #   ghtMode()
    
    # FUNCTION:
    
    # Find Maximum: 
    min = qght(0.01, beta, delta, mu, nu)
    max = qght(0.99, beta, delta, mu, nu)
    ans = optimize(f = dght, interval = c(min, max), 
        beta = beta, delta = delta, mu = mu, nu = nu, 
        maximum = TRUE, tol = .Machine$double.eps)$maximum
        
    # Return Value:
    ans
} 


################################################################################

