
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
# FUNCTION:               DESCRIPTION:
#  gridVector              Creates from two vectors rectangular grid coordinates
################################################################################


gridVector <-  
function(x, y = NULL)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Creates from two vectors rectangular grid coordinates
    
    # Arguments:
    #   x, y - two numeric vectors of length m and n which span the 
    #   rectangular grid of size m times n.
    
    # Details: 
    #   The two vectors x and y span a rectangular grid with nx=length(x) 
    #   times ny=length(y) points which are returned as a matrix of size
    #   (nx*ny) times 2.
    
    # Value:
    #   returns a list with two elements X and Y each of length m 
    #   times n
    
    # Example:
    #   > gridVector(1:3, 1:2)
    #             [,1] [,2]
    #       [1,]    1    1
    #       [2,]    2    1
    #       [3,]    3    1
    #       [4,]    1    2
    #       [5,]    2    2
    #       [6,]    3    2

    # FUNCTION: 
    
    # DW:
    if (is.null(y)) y = x
    
    # Prepare for Input:
    nx = length(x)
    ny = length(y)
    xoy = cbind(rep(x, ny), as.vector(matrix(y, nx, ny, byrow = TRUE)))
    X = matrix(xoy, nx * ny, 2, byrow = FALSE)
    
    # Return Value:
    list(X = X[,1], Y = X[,2])
}


################################################################################

