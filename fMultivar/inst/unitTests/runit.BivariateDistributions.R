
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
#  grid2d                Returns from two vectors x-y grid coordinates
#  density2d             Returns 2D Kernel Density Estimates
#  hist2d                Returns 2D Histogram Counts
# FUNCTION:             BIVARIATE DISTRIBUTIONS:
#  pnorm2d               Computes bivariate Normal probability function
#  dnorm2d               Computes bivariate Normal density function
#  rnorm2d               Generates bivariate normal random deviates
#  pcauchy2d             Computes bivariate Cauchy probability function
#  dcauchy2d             Computes bivariate Cauchy density function
#  rcauchy2d             Generates bivariate Cauchy random deviates
#  pt2d                  Computes bivariate Student-t probability function
#  dt2d                  Computes bivariate Student-t density function
#  rt2d                  Generates bivariate Student-t random deviates
# FUNCTION:             ELLIPTICAL DISTRIBUTIONS:
#  delliptical2d         Computes density for elliptical distributions
# REQUIREMENTS:
#  fBasics::.perspPlot
#  fBasics::.contourPlot
################################################################################


test.grid2d =
function()
{
    # Grid Data:
    grid2d(x = (0:10)/10)
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.density2d =
    function()
{
    # Data:
    z <- rnorm2d(1000)
    
    # Density:
    D = density2d(x = z[, 1], y = z[, 2])
    .perspPlot(D)
    .contourPlot(D)
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.hist2d =
    function()
{
    # Data:
    z <- rnorm2d(1000)
   
    # Histogram:
    H <- hist2d(x = z[, 1], y = z[, 2])
    .perspPlot(H)
    .contourPlot(H)

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.norm2d =
    function()
{
    #  pnorm2d - Computes bivariate Normal probability function
    #  dnorm2d - Computes bivariate Normal density function
    #  rnorm2d - Generates bivariate normal random deviates
    
    # Normal Density:
    x = (-40:40)/10
    X = grid2d(x)
    z = dnorm2d(X$x, X$y)
    Z = list(x = x, y = x, z = matrix(z, ncol = length(x)))
    .perspPlot(Z)
    .contourPlot(Z)
    
    # Normal Density, rho = 0.5:
    x = (-40:40)/10
    X = grid2d(x)
    z = dnorm2d(X$x, X$y, rho = 0.5)
    Z = list(x = x, y = x, z = matrix(z, ncol = length(x)))
    .perspPlot(Z)
    .contourPlot(Z)
 
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.cauchy2d =
    function()
{
    #  pcauchy2d - Computes bivariate Cauchy probability function
    #  dcauchy2d - Computes bivariate Cauchy density function
    #  rcauchy2d - Generates bivariate Cauchy random deviates
    
    # Cauchy Density:
    x = (-40:40)/10
    X = grid2d(x)
    z = dcauchy2d(X$x, X$y)
    Z = list(x = x, y = x, z = matrix(z, ncol = length(x)))
    .perspPlot(Z)
    .contourPlot(Z)
    
    # Cauchy Density, rho = 0.5:
    x = (-40:40)/10
    X = grid2d(x)
    z = dcauchy2d(X$x, X$y, rho = 0.5)
    Z = list(x = x, y = x, z = matrix(z, ncol = length(x)))
    .perspPlot(Z)
    .contourPlot(Z)
 
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.t2d =
    function()
{
    #  pt2d - Computes bivariate Student-t probability function
    #  dt2d - Computes bivariate Student-t density function
    #  rt2d - Generates bivariate Student-t random deviates
    
    # Student Density:
    x = (-40:40)/10
    X = grid2d(x)
    z = dt2d(X$x, X$y, nu = 4)
    Z = list(x = x, y = x, z = matrix(z, ncol = length(x)))
    .perspPlot(Z)
    .contourPlot(Z)
    
    # Student Density, rho = 0.5:
    x = (-40:40)/10
    X = grid2d(x)
    z = dt2d(X$x, X$y, rho = 0.5, nu = 4)
    Z = list(x = x, y = x, z = matrix(z, ncol = length(x)))
    .perspPlot(Z)
    .contourPlot(Z)
 
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.delliptical2d =
    function()
{ 
    # Settings:
    xy = grid2d((-50:50)/10)
    
    # Contour Plots:
    par(mfrow = c(1, 1))
    contour(delliptical2d(xy, rho = 0.75, param = NULL, 
        type = "norm", output = "list"), main = "norm")
    contour(delliptical2d(xy, rho = 0.75, param = NULL, 
        type = "cauchy", output = "list"), main = "cauchy")
    contour(delliptical2d(xy, rho = 0.75, param = 4, 
        type = "t", output = "list"), main = "t")
    contour(delliptical2d(xy, rho = 0.75, param = NULL, 
        type = "laplace", output = "list"), main = "laplace")
    contour(delliptical2d(xy, rho = 0.75, param = NULL, 
        type = "kotz", output = "list"), main = "kotz")
    contour(delliptical2d(xy, rho = 0.75, param = NULL, 
        type = "epower", output = "list"), main = "epower")
        
    # Return Value:
    return()    
}


################################################################################
    
