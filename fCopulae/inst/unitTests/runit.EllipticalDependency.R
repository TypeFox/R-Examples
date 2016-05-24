
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
# FUNCTION:                  ELLIPTICAL COPULAE DEPENDENCE MASURES:
#  ellipticalTau              Computes Kendall's tau for elliptical copulae
#  ellipticalRho              Computes Spearman's rho for elliptical copulae
# FUNCTION:                  ELLIPTICAL COPULAE TAIL COEFFICIENT:
#  ellipticalTailCoeff        Computes tail dependence for elliptical copulae
#  ellipticalTailPlot         Plots tail dependence function
################################################################################


test.ellipticalTau = 
function()
{
    # Computes Kendall's tau for elliptical copulae
    args(ellipticalTau)
    
    ellipticalTau(rho = 0.5)
    ellipticalTau(rho = c(-0.5, 0, 0.5))
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.ellipticalRho = 
function()
{    
    # Computes Spearman's rho for elliptical copulae
    args(ellipticalRho)
    
    ellipticalRho(0.5)
    ellipticalRho(rho = c(-0.5, 0, 0.5))
 
    # Return Value:
    return()    
}


################################################################################


test.ellipticalTailCoeff = 
function()
{
    # Lower - Upper ----
    
    # Tail Coefficient - Using Default Parameters:
    Type = c("norm", "cauchy", "t") 
    for (type in Type) {
        ans = ellipticalTailCoeff(rho = 0.5, type = type)   
        print(ans)
        cat("\n")
    }
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.ellipticalTailPlot = 
function()
{    
    # Arguments:
    # ellipticalTailPlot(param = NULL, type = c("norm", "cauchy", "t"), 
    #   tail = c("Lower", "Upper")) 
    
    # Plot - Be patient, plotting takes some time ...
    Type = c("norm", "cauchy", "t") 
    for (type in Type) {
        par(mfrow = c(2, 2), cex = 0.7)
        ellipticalTailPlot(type = type)  
        ellipticalTailPlot(type = type, tail = "Lower")   
    }           
    
    # Return Value:
    return()    
}

  
################################################################################
   
