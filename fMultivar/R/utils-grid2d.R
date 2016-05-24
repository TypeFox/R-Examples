
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
#  grid2d                Returns from two vectors x-y grid coordinates
################################################################################


grid2d <- 
    function(x = (0:10)/10, y = x)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Creates from two vectors x-y grid coordinates
    
    # Arguments:
    #   x, y - two numeric vectors defining the x and y coordinates.
    
    # Value:
    #   returns a list with two vectors named $x and $y spanning the 
    #   grid defined by the coordinates x and y.
    
    # Example:
    #   > grid2d(1:3, 1:2)
    #       $x
    #       [1] 1 2 3 1 2 3
    #       $y
    #       [1] 1 1 1 2 2 2

    # FUNCTION: 
    
    # Prepare for Input:
    nx  = length(x)
    ny  = length(y)
    xoy = cbind(rep(x, ny), as.vector(matrix(y, nx, ny, byrow = TRUE)))
    XY  = matrix(xoy, nx * ny, 2, byrow = FALSE)
    
    # Return Value:
    list(x = XY[, 1], y = XY[, 2])
}



################################################################################

