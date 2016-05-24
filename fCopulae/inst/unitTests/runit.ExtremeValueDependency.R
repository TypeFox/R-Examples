
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
#  evTau                      Returns Kendall's tau for extreme value copulae
#  evRho                      Returns Spearman's rho for extreme value copulae
# FUNCTION:                  EXTREME VALUE COPULAE TAIL DEPENDENCE:
#  evTailCoeff                Computes tail dependence for extreme value copulae
#  evTailCoeffSlider          Plots extreme value tail dependence function
#################################################################################


test.evTau = 
function()
{
    # Arguments:
    #   evTau(param = NULL, type = evList(), alternative = FALSE)
    
    # Tau:
    for (type in evList()) {
        ans = evTau(type = type)
        cat("\n")
        print(ans)
    }
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.evRho = 
function()
{
    # Arguments:
    #   evRho(param = NULL, type = evList(), alternative = FALSE) 
    
    # Rho:
    for (type in evList()) {
        ans = evRho(type = type)
        cat("\n")
        print(type)
        print(ans)
    }
    
    # Return Value:
    return()    
}

   
################################################################################
    

test.evTailCoeff = 
function()
{
    # Arguments:
    #   evTailCoeff(param = NULL, type = evList())
    
    # Tail Coefficient:
    for (type in evList()) {
        ans = evTailCoeff(type = type)
        cat("\n")
        print(type)
        print(ans)
    }
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------
    

test.evTailCoeffSlider = 
function()
{
    # Arguments:
    #   evTailCoeffSlider(B = 10) 
    
    # Try Slider:
    # evTailCoeffSlider()                                                   
    NA
    
    # Return Value:
    return()    
}

  
################################################################################

