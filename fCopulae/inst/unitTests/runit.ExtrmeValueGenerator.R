
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
# FUNCTION:                 EXTREME VALUE COPULAE PARAMETER:
#  evList                    Returns list of implemented extreme value copulae
#  evParam                   Sets Default parameters for an extreme value copula
#  evRange                   Returns the range of valid parameter values
#  evCheck                   Checks if parameters are in the valid range
# FUNCTION:                 EXTREME VALUE COPULAE GENERATOR FUNCTION:
#  Afunc                     Computes Dependence function
#  AfuncSlider               Displays interactively dependence function
#################################################################################


test.evList = 
function()
{
    # Arguments:
    #   evList()
    
    # List:
    evList()
    # c("gumbel", "galambos", "husler.reiss", "tawn", "bb5")

    # Return Value:
    return()    
}

    
# ------------------------------------------------------------------------------


test.evParam = 
function()
{
    # Arguments:
    #   evParam(type = evList())
    
    # Parameters:
    for (type in evList()) {
        cat("\n")
        print(unlist(evParam(type)))
    }

    # Return Value:
    return()    
}

    
# ------------------------------------------------------------------------------


test.evRange = 
function()
{
    # Arguments:
    # evRange(type = evList()) 

    # Range:
    for (type in evList()) {
        cat("\n")
        print(evRange(type))
    }

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.evCheck = 
function()
{
    # Arguments:
    # evCheck(type = evList()) 
    
    # Check:
    for (type in evList()) {
        cat("\n")
        param = evParam(type)$param
        print(evCheck(param))
    }
    
    # Return Value:
    return()    
}

   
################################################################################


test.Afunc = 
function()
{
    # Arguments:
    #   Afunc(x, param = NULL, type = evList()
    
    # Afunc:
    x = (0:10)/10
    for (type in evList()) {
        cat("\n")
        print(type)
        print(Afunc(x, type = type))
    }
    
    # Return Value:
    return()    
}

   
# ------------------------------------------------------------------------------

    
test.AfuncSlider = 
function()
{
    # Arguments:
    #   AfuncSlider()
    
    # Try Slider:
    # AfuncSlider()
    NA
    
    # Return Value:
    return()    
}

  
################################################################################

