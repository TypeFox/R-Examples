
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
# FUNCTION:                  EXTREME VALUE COPULA PARAMETER FITTING:
#  evCopulaSim                Simulates bivariate extreme value copula
#  evCopulaFit                Fits the paramter of an extreme value copula
#################################################################################


test.evCopulaSim = 
function()
{
    # Arguments:
    # evCopulaSim(n, param = NULL, type = evList()) 

    # Simulate Random Variates:
    for (type in evList()) {
        ans = evCopulaSim(5, type = type)
        cat("\n")
        print(type)
        print(ans)
    }
    
    # Return Value:
    return()    
}

    
# ------------------------------------------------------------------------------


test.evCopulaFit = 
function()
{
    # Arguments:
    #   evCopulaFit(u, v = NULL, type = evList(), ...) 

    # Random Variates:
    set.seed(4711)
    type = "gumbel"
    R = evCopulaSim(500, param = NULL, type = type)
    Index =  which(is.na(R[,2]))
    R = R[-Index, ] 
    
    # Fit:
    ### evCopulaFit(u = R, type = type)                                  # Check
    
    # Fit:
    ### evCopulaFit(u = R[, 1], v = R[, 2], type = type)                 # Check
    
    # Return Value:
    return()    
}

  
################################################################################

