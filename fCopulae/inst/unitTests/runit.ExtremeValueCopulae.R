
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


################################################################################
# FUNCTION:                  EXTREME VALUE COPULAE RANDOM VARIATES:
#  revCopula                  Generates extreme value copula random variates 
#  revSlider                  isplays interactively plots of random variates
# FUNCTION:                  EXTREME VALUE COPULAE PROBABILIY:
#  pevCopula                  Computes extreme value copula probability
#  pevSlider                  Displays interactively plots of probability
# FUNCTION:                  EXTREME VALUE COPULAE DENSITY:
#  devCopula                  Computes extreme value copula density
#  devSlider                  Displays interactively plots of density
################################################################################


test.revCopula = 
function()
{
    # Arguments:
    # revCopula(n, param = NULL, type = evList())

    # Random Variates - Check all Types:
    for (type in evList()) {
        R = revCopula(n = 5, param = NULL, type = type)
        cat("\n")
        print(type)
        print(R)
    }
     
    # Tawn Copula:
    revCopula(n = 5, param = NULL, type = "tawn")
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.revSlider = 
function()
{
    # Arguments:
    # revSlider(B = 10)
    
    # Try Slider()
    # revSlider()                         # CHECK !!!
    NA
    
    # Return Value:
    return()    
}

    
################################################################################


test.pevCopula = 
function()
{
    # Arguments:
    # pevCopula(u = 0.5, v = u, param = NULL, type = evList(), 
    #   output = c("vector", "list"), alternative = FALSE) 

    # Random Variates - Check all Types:
    for (type in evList()) {
        R = pevCopula(u = grid2d(), param = NULL, type = type, output = "list")
        cat("\n")
        print(type)
        print(R)
    }
     
    # Tawn Copula:
    revCopula(n = 5, param = NULL, type = "tawn")
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.pevSlider = 
function()
{
    # Arguments:
    # pevSlider(type = c("persp", "contour"), B = 10)

    # Try Perspective Slider:
    # pevSlider("persp")
    NA
    
    # Try Contour Slider:
    # pevSlider("contour")
    NA
    
    # Return Value:
    return()    
}


################################################################################


test.devCopula = 
function()
{
    # Arguments:
    # devCopula(u = 0.5, v = u, param = NULL, type = evList(), 
    #   output = c("vector", "list"), alternative = FALSE)

    # Random Variates - Check all Types:
    for (type in evList()) {
        R = devCopula(u = grid2d(), param = NULL, type = type, output = "list")
        cat("\n")
        print(type)
        print(R)
    }                                           # CHECK Border !!!!
     
    # Tawn Copula:
    revCopula(n = 5, param = NULL, type = "tawn")
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.devSlider = 
function()
{
    # Arguments:
    # devSlider(type = c("persp", "contour"), B = 10)

    # Try Perspective Slider:
    # devSlider("persp")
    NA 
    
    # Try Contour Slider:
    # devSlider("contour")
    NA
    
    # Return Value:
    return()    
}

  
################################################################################

