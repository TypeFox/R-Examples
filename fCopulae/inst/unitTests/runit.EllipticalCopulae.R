
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
# FUNCTION:                  ELLIPTICAL COPULAE RANDOM DEVIATES:
#  rellipticalCopula          Generates elliptical copula variates
#  rellipticalSlider          Interactive plots of random variates
# FUNCTION:                  ELLIPTICAL COPULAE PROBABILITY:
#  pellipticalCopula          Computes elliptical copula probability
#  pellipticalSlider          Interactive plots of probability
# FUNCTION:                  ELLIPTICAL COPULAE DENSITY:
#  dellipticalCopula          Computes elliptical copula density 
#  dellipticalSlider          Interactive plots of density
################################################################################


test.rellipticalCopula = 
    function()
{
    # Random Number Generator:
    R <- rellipticalCopula(1000, type = "norm")
    plot(R, pch = 19, col = "steelblue", main = "norm")
    grid()
    
    R <- rellipticalCopula(1000, type = "cauchy")
    plot(R, pch = 19, col = "steelblue", main = "cauchy")
    grid()
    
    R <- rellipticalCopula(1000, type = "t") 
    plot(R, pch = 19, col = "steelblue", main = "t-default")
    grid()
    
    R <- rellipticalCopula(1000, param = c(nu = 3), type = "t")
    plot(R, pch = 19, col = "steelblue", main = "t3")
    grid()
    
    R <- rellipticalCopula(1000, param = 3, type = "t")
    plot(R, pch = 19, col = "steelblue", main = "t3")
    grid()
    
    # The remaining copulae are not yet implemented ...
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.rellipticalSlider <-  
    function()
{   
    # Try Slider:
    # rellipticalSlider()
    NA
    
    # Return Value:
    return()    
}


################################################################################


test.pellipticalCopula <- 
    function()
{ 
    # Arguments ?
    # pellipticalCopula(u = 0.5, v = u, rho = 0.75, param = NULL, 
    #   type = ellipticalList(), output = c("vector", "list"), border = TRUE) 
    
    # Use Default Settings:
    par (mfrow = c(1, 1))
    for (type in ellipticalList()) {
        UV <- grid2d()
        p <- pellipticalCopula(u = UV, rho = 0.75, type = type, output = "list")   
        print(type)
        persp(p, main = type, theta = -40, phi = 30, col = "steelblue", 
            ps = 9, xlab = "u", ylab = "v", zlab = "C")
    }
   
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.pellipticalSlider = 
function()
{       
    # Arguments:
    # pellipticalSlider(type = c("persp", "contour"), B = 20) 
 
    # Try Perspective Slider:
    # pellipticalSlider()
    NA
    
    # Try Contour Slider:
    # pellipticalSlider("contour")
    NA
   
    # Return Value:
    return()    
}


################################################################################


test.dellipticalCopula = 
function()
{  
    # Arguments ?
    # dellipticalCopula(u = 0.5, v = u, rho = 0.75, param = NULL, 
    #   type = ellipticalList(), output = c("vector", "list"), border = TRUE) 

    # Use Default Settings:
    par (mfrow = c(1, 1))
    for (type in ellipticalList()) {
        UV = grid2d()
        d = dellipticalCopula(u = UV, rho = 0.75, type = type, output = "list")   
        print(type)
        persp(d, main = type, theta = -40, phi = 30, col = "steelblue", 
            ps = 9, xlab = "u", ylab = "v", zlab = "c")
    }
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.dellipticalSlider = 
function()
{  
    # Arguments:
    # dellipticalSlider(type = c("persp", "contour"), B = 20) 
    
    # Try Perspective Slider:
    # dellipticalSlider()
    NA
    
    # Try Contour Slider:
    # dellipticalSlider("contour")
    NA
    
    # Return Value:
    return()    
}

  
################################################################################
   
