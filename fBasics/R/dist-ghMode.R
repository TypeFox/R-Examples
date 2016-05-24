
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
#  ghMode                Computes mode of the generalized hyperbolic DF
################################################################################


ghMode <- 
function(alpha = 1, beta = 0, delta = 1, mu = 0, lambda = -1/2)
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Computes the mode of the Generalized Hyperbolic PDF
  
    # FUNCTION:
    
    # Find Maximum: 
    min = qgh(0.01, alpha, beta, delta, mu, lambda)
    max = qgh(0.99, alpha, beta, delta, mu, lambda)
    ans = optimize(f = dgh, interval = c(min, max), 
        alpha = alpha, beta = beta, delta = delta, mu = mu, lambda = lambda, 
        maximum = TRUE, tol = .Machine$double.eps)$maximum
        
    # Return Value:
    ans
} 


################################################################################

