
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
# FUNCTION                   KENDALL'S TAU AND SPEARMAN'S RHO:
#  archmTau                   Returns Kendall's tau for Archemedean copulae
#  archmRho                   Returns Spearman's rho for Archemedean copulae
# FUNCTION:                  ARCHIMEDEAN COPULAE TAIL COEFFICIENT:
#  archmTailCoeff             Computes tail dependence for Archimedean copulae
#  archmTailPlot              Plots Archimedean tail dependence function
################################################################################


test.archmTau = 
function()
{
    # Arguments:
    # archmTau(alpha = NULL, type = archmList(), lower = 1e-10)
    
    # Tau:
    for (type in archmList()) {
        ans = archmTau(type = type)
        cat("\n")
        print(type)
        print(ans)
    }
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.archmRho = 
function()
{
    # Arguments:
    # archmRho(alpha = NULL, type = archmList(), 
    #   method = c("integrate2d", "adapt"), error = 1e-05)
    
    # Rho:
    for (type in archmList()) {
        ans = archmRho(alpha = NULL, type = type, 
            method = "integrate2d", error = 1e-5)
        cat("\n")
        print(type)
        print(ans)
    }
    
    # Return Value:
    return()    
}


################################################################################


test.archmTailCoeff = 
function()
{
    # Arguments:
    # archmTailCoeff(alpha = NULL, type = archmList())
    
    # Tail Coefficient:
    for (type in archmList()) {
        ans = archmTailCoeff(alpha = NULL, type = type)
        cat("\n")
        print(type)
        print(ans)
    }
        
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.archmTailPlot = 
function()
{
    # Arguments:
    # archmTailPlot(alpha = NULL, type = archmList(), 
    #   tail = c("Upper", "Lower")) 
    
    # Lower Tail Coefficient Plot:
    par(mfrow = c(2, 2), cex = 0.7)
    for (type in archmList()) {
        print(type)
        archmTailPlot(alpha = NULL, type = type, tail = "Upper")
    }
    
    # Upper Tail Coefficient Plot:
    for (type in archmList()) {
        print(type)
        archmTailPlot(alpha = NULL, type = type, tail = "Lower")
    }
    
    # Return Value:
    return()    
}
 
  
################################################################################

