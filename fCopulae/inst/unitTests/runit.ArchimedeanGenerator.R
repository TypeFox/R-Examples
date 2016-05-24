
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
# FUNCTION:                 ARCHIMEDEAN COPULAE PARAMETER:
#  evList                    Returns list of implemented extreme value copulae
#  archmParam                Sets Default parameters for an extreme value copula
#  archmRange                Returns the range of valid alpha values
#  archmCheck                Checks if alpha is in the valid range
# FUNCTION:                 ARCHIMEDEAN COPULAE PHI GENERATOR:
#  Phi                       Computes Archimedean Phi, inverse and derivatives
#  PhiSlider                 Displays interactively generator function
# FUNCTION:                 ARCHIMEDEAN DENSITY K GENERATOR:
#  Kfunc                     Computes Archimedean Density Kc and its Inverse
#  KfuncSlider               Displays interactively the density and concordance
################################################################################


test.archmList = 
function()
{

    # Arguments:
    # archmList()
    
    # List:
    archmList()
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.archmParam = 
function()
{
    # Arguments:
    # archmParam(type = archmList()) 
    
    # Parameters:
    for (type in archmList()) {
        cat("\n")
        print(unlist(archmParam(type)))
    }
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.archmRange = 
function()
{
    # Arguments:
    # archmRange(type = archmList(), B = Inf)
    
    # Range:
    for (type in archmList()) {
        cat("\n")
        print(archmRange(type))
    }

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.archmCheck = 
function()
{
    # Arguments ?
    # archmCheck(alpha, type = archmList())
    
    # Check:
    for (type in archmList()) {
        cat("\n")
        print(archmCheck(archmParam(type)$param))
    }

    # Return Value:
    return()    
}


################################################################################


test.Phi = 
function()
{
    # Arguments:
    # Phi(x, alpha = NULL, type = archmList(), inv = FALSE, deriv = paste(0:2))

    # Call Generator Function Phi:
    for (type in paste(1:22)) {
        print(Phi(x = 0.5, type = type, inv = TRUE, deriv = "0"))
        cat("\n")
    }
    for (type in paste(1:22)) {
        print(Phi(x = 0.5, type = type, inv = TRUE, deriv = "1"))
        cat("\n")
    }
    for (type in paste(1:22)) {
        print(Phi(x = 0.5, type = type, inv = TRUE, deriv = "2"))
        cat("\n")
    }
    
    for (type in paste(1:22)) {
        print(Phi(x = 0.5, type = type, inv = FALSE, deriv = "0")) 
        cat("\n")
    }
    for (type in paste(1:22)) {
        print(Phi(x = 0.5, type = type, inv = FALSE, deriv = "1"))
        cat("\n")
    }
    for (type in paste(1:22)) {
        print(Phi(x = 0.5, type = type, inv = FALSE, deriv = "2"))
        cat("\n")
    }
   
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.PhiSlider = 
function()
{
    # Arguments:
    # PhiSlider()
    
    # Try Slider:
    # PhiSlider()
    NA
    
    # Return Value:
    return()    
}


################################################################################


test.Kfunc = 
function()
{
    # Arguments:
    # Kfunc(x, alpha = NULL, type = archmList(), inv = FALSE, lower = 1e-08)
    
    # Call Generator Function Phi:
    for (type in paste(1:22)) {
        print(Kfunc(x = 0.5, inv = FALSE))
        cat("\n")
    }
    for (type in paste(1:22)) {
        print(Kfunc(x = 0.5, inv = TRUE))
        cat("\n")
    }
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------
 
    
test.KfuncSlider = 
function()
{
    # Arguments:
    # KfuncSlider()
    
    # Try Slider:
    # KfuncSlider()
    NA
    
    # Return Value:
    return()    
}
   

################################################################################

