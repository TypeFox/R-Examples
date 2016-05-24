
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
# FUNCTION:                 BIVARIATE GRIDDED INTERPOLATION:
#  linearInterp             Interpolates linearly irregularly spaced data points
#  linearInterpp            Interpolates linearly pointwise
################################################################################


linearInterp <-
  function(x, y = NULL, z = NULL, gridPoints = 21,
           xo = seq(min(x), max(x), length = gridPoints),
           yo = seq(min(y), max(y), length = gridPoints))
  {
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Interpolates Linearly Irregularly Distributed Data Points
    
    # Arguments:
    #   x, y, z - either three numeric vectors of equal length or if
    #       y and z are NULL, a list with entries x, y, a, or named
    #       data.frame with x in the first, y in the second, and z in
    #       the third column.
    #   gridPoints - number of grid points in x and y direction.
    #   xo, yo, a sequence of data points spanning the grid
    
    # Note:
    #   Extrapolation is not possible in the case of linear interpolation.
    
    # Value:
    #   A list with three elements, $x and $y which are vectors of length
    #   'gridPoints' and $z which is a matrix of size 'gridPoints^2'.
    
    # Requirements:
    #   akima Builtin Fortran Code.
    
    # Example:
    #   set.seed(1953)
    #   x = runif(999)-0.5; y = runif(999)-0.5; z = cos(2*pi*(x^2+y^2))
    #   ans = linearInterp(x, y, z)
    #   persp(ans, theta = -50, phi = 30, col = "steelblue")
    
    # Note:
    #   Uses Fortran akima Builtin
    
    # FUNCTION:
    
    if (!require(akima, quietly = TRUE))
      stop("\n -- Package akima not available -- \n\n")
    
    # Arguments:
    if (is.list(x)) x = matrix(unlist(x), ncol = 3)
    if (is.data.frame(x)) x = as.matrix.data.frame(x)
    if (is.matrix(x)) {
      z = x[, 3]
      y = x[, 2]
      x = x[, 1]
    }
    
    # Interpolation:
    ans <- akima::interp.old(x, y, z, xo, yo, ncp = 0, extrap = FALSE,
                             duplicate = "median", dupfun = NULL)
    colnames(ans$z) <- as.character(signif(ans$x, round(log(gridPoints), 0)))
    rownames(ans$z) <- as.character(signif(ans$y, round(log(gridPoints), 0)))
    class(ans) <- "gridData"
    
    # Return Value:
    ans
  }


# ------------------------------------------------------------------------------


linearInterpp <-
  function(x, y = NULL, z = NULL, xo, yo)
  {
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Interpolates Linearly Irregularly Distributed Data Points
    
    # Arguments:
    #   x, y, z - either three numeric vectors of equal length or if
    #       y and z are NULL, a list with entries x, y, a, or named
    #       data.frame with x in the first, y in the second, and z in
    #       the third column.
    #   gridPoints - number of grid points in x and y direction.
    #   xo, yo, a sequence of data points for pointwise interpolation
    
    # Note:
    #   Extrapolation is not possible in the case of linear interpolation.
    
    # Value:
    #   A list with three elements, $x and $y which are vectors of length
    #   'gridPoints' and $z which is a matrix of size 'gridPoints^2'.
    
    # Requirements:
    #   akima Builtin Fortran Code.
    
    # Example:
    #   set.seed(1953)
    #   x = runif(999)-0.5; y = runif(999)-0.5; z = cos(2*pi*(x^2+y^2))
    #   ans = linearInterpp(x, y, z, c(mean(x), 0, 100), c(mean(y), 0, 100))
    #   persp(ans, theta = -50, phi = 30, col = "steelblue")
    
    # Note:
    #   Uses Fortran akima Builtin
    
    # FUNCTION:
    
    if (!require(akima, quietly = TRUE))
      stop("\n -- Package akima not available -- \n\n")
    
    # Arguments:
    if (is.list(x)) x = matrix(unlist(x), ncol = 3)
    if (is.data.frame(x)) x = as.matrix.data.frame(x)
    if (is.matrix(x)) {
      z = x[, 3]
      y = x[, 2]
      x = x[, 1]
    }
    
    # Interpolation:
    interpp.old <- eval(parse(text=paste0("akima",":::","interpp.old")))
    ans <- interpp.old(x, y, z, xo, yo, ncp = 0, extrap = FALSE,
                       duplicate = "median", dupfun = NULL)
    ans <- data.frame(x = ans$x, y = ans$y, z = ans$z)
    
    # Return Value:
    ans
  }


################################################################################
