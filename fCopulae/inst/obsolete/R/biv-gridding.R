
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
# You should have received A copy of the GNU Library General 
# Public License along with this library; if not, write to the 
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
# MA  02111-1307  USA

# Copyrights (C)
# for this R-port: 
#   1999 - 2007, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


# fEcofin::4A-BivariateGridding.R
################################################################################
# FUNCTION:                 GRID DATA:
#  gridData                  Generates grid data set
#  persp.gridData            Generates perspective plot from a grid data object
#  contour.gridData          Generates contour plot from a grid data object
################################################################################


################################################################################
# FUNCTION:                 GRID DATA:
#  gridData                  Generates grid data set
#  persp.gridData            Generates perspective plot from a grid data object
#  contour.gridData          Generates contour plot from a grid data object


gridData =
function(x = (-10:10)/10, y = x, z = outer(x, y, function(x, y) (x^2+y^2)) )
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Generates a grid data set
    
    # Arguments:
    #   x, y -  two numeric vectors of grid pounts
    #   z -  a numeric matrix or any other rectangular object which can 
    #       be transformed by the function 'as.matrix' into a matrix
    #       object.
    
    # Example:
    #   persp(as.gridData())
    
    # FUNCTION:
    
    # Grid Data:
    data = list(x = x, y = y, z = as.matrix(z))
    class(data) = "gridData"
    
    # Return Value:
    data
}
 

# ------------------------------------------------------------------------------


persp.gridData =
function(x, theta = -40, phi = 30, col = "steelblue", ticktype = "detailed",
...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   S3 method to generate a perspective plot from a grid data object
    
    # Example:
    #   x = y = seq(-10, 10, length = 30)
    #   z = outer(x, y, function(x, y) { r <- sqrt(x^2+y^2); 10 * sin(r)/r })
    #   data = list(x = x, y = y, z = z)
    #   class(data) = "gridData"
    #   persp(data)
    
    # FUNCTION:
    
    # Grid Data:
    class(x) = "default"
    persp(x, theta = theta, phi = phi, col = col, ticktype = ticktype, ...) 
        
    # Return Value:
    invisible(NULL)    
}


# ------------------------------------------------------------------------------


contour.gridData =
function(x, addImage = TRUE, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   S3 method to generate a contour plot from a grid data object
    
    # Example:
    #   x = y = seq(-10, 10, length = 30)
    #   z = outer(x, y, function(x, y) { r <- sqrt(x^2+y^2); 10 * sin(r)/r })
    #   data = list(x = x, y = y, z = z)
    #   class(data) = "gridData"
    #   contour(data)

    # FUNCTION:

    # Grid Data:
    class(x) = "default"
    if (addImage) image(x, ...)
    contour(x, add = addImage, ...)
    box()   
    
    # Return Value:
    invisible(NULL)
}


################################################################################

