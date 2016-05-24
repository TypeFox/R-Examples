
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


################################################################################
# FUNCTION:             DESCRIPTION:
#  integrate2d           Integrates over a two dimensional unit square
################################################################################


integrate2d <-  
    function(fun, error = 1.0e-5, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   2-dimension quadrature rule on [0,1]^2 
    
    # Arguments:
    #   fun - function to be integrated. The first argument requests
    #       the x values, the second the y values, and the remaining
    #       are reserved for additional parameters.
    #   ... - parameters passed to the function to be integrated
    
    # Details:
    #   see: Abramowitz and Stegun, p. 892
    
    # FUNCTION:
    
    # Estimate a reasonable number of subintervals:
    H <- sqrt(sqrt(error))
    n <- ceiling(1/H + 1)
    blocks <- ceiling(log(n+1)/log(2))
    n <- 2^blocks-1
    h <- 1/(n-1)
    
    # The error will be of order h^4:
    error <- h^4
    
    # Create all grid coordinates:
    x <- y <- h*seq(1, n-1, by = 2)
    nx <- ny <- length(x)
    xoy <- cbind(rep(x, ny), as.vector(matrix(y, nx, ny, byrow = TRUE)))
    XY <- matrix(xoy, nx * ny, 2, byrow = FALSE)
    
    # The integration rule:      
    rule <- function(x, h, ...) {         
        X = x[1] + h*c(  0, -1, -1,  1,  1, -1,  1,  0,  0)
        Y = x[2] + h*c(  0, -1,  1, -1,  1,  0,  0, -1,  1)
        W =  c( 16,  1,  1,  1,  1,  4,  4,  4,  4)/36
        ans = sum( W * fun(X, Y, ...) )
    }
     
    # Result:
    ans <- (4*h^2)*sum(apply(XY, 1, rule, h = h, ...))
    
    # Return Value:
    list(value = ans, error = error, points = n)
}


################################################################################

