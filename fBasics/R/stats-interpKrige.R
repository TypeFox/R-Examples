
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


################################################################################
# FUNCTION:                DESCRIPTION:
#  krigeInterp              Kriges irregularly distributed data points
# REQUIRE:                 DESCRIPTION:
#  spatial                  Functions for Kriging and Point Pattern 
################################################################################


krigeInterp <- 
  function(x, y = NULL, z = NULL, gridPoints = 21,
           xo = seq(min(x), max(x), length = gridPoints),
           yo = seq(min(y), max(y), length = gridPoints), extrap = FALSE,
           polDegree = 6)
  {
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Kriges Irregularly Distributed Data Points
    
    # Arguments:
    #   x, y, z - either three numeric vectors of equal length or if
    #       y and z are NULL, a list with entries x, y, a, or named
    #       data.frame with x in the first, y in the second, and z in
    #       the third column.
    #   gridPoints - number of grid points in x and y direction.
    #   xo, yo, a sequence of data points spanning the grid
    #   extrap - a logical, if TRUE then the data points are extrapolated.
    #   polDegree - polynomial degree, an integer ranging between 1 and 6.
    
    # Value:
    #   A list with three elements, $x and $y which are vectors of length
    #   'gridPoints' and $z which is a matrix of size 'gridPoints^2'.
    
    # Example:
    #   x <- runif(999)-0.5; y = runif(999)-0.5; z = cos(2*pi*(x^2+y^2))
    #   require(spatial)
    #   ans <- krigeInterp(x, y, z, extrap = FALSE)
    #   persp(ans, theta = -50, phi = 30, col = "steelblue")
    
    # Note:
    #   Requires recommended R Package "spatial"
    
    # FUNCTION:
    
    if (!require(spatial, quietly = TRUE))
      stop("\n -- Package spatial not available -- \n\n")
    
    # Arguments:
    if (is.list(x)) x <- matrix(unlist(x), ncol = 3)
    if (is.data.frame(x)) x <- as.matrix.data.frame(x)
    if (is.matrix(x)) {
      z = x[, 3]
      y = x[, 2]
      x = x[, 1]
    }
    
    # Interpolate:
    krige <- spatial::surf.gls(np = polDegree, covmod = spatial::expcov,
                               x = x, y = y, z = z, d = 0.5, alpha = 1)
    ans <- spatial::prmat(krige,
                          xl = min(xo), xu = max(xo), yl = min(yo), yu = max(yo),
                          n = gridPoints-1)
    
    # Extrapolate ?
    # - this should be done more efficiently
    if (!extrap) {
      E <- akimaInterp(x = x, y = y, z = z, gridPoints = gridPoints,
                       extrap = extrap)
      ans$z[is.na(E$z)] = NA
    }
    class(ans) <- "gridData"
    
    # Return Value:
    ans
  }


################################################################################

