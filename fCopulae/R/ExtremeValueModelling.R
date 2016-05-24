
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
# FUNCTION:                  EXTREME VALUE COPULA PARAMETER FITTING:
#  evCopulaSim                Simulates bivariate extreme value copula
#  evCopulaFit                Fits the paramter of an extreme value copula
################################################################################


################################################################################
# FUNCTION:                  EXTREME VALUE COPULA PARAMETER FITTING:
#  evCopulaSim                Simulates bivariate extreme value copula
#  evCopulaFit                Fits the paramter of an extreme value copula


evCopulaSim = 
function(n, param = NULL, type = evList())
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Simulates bivariate extreme value Copula
    
    # FUNCTION:
    
    # Match Arguments:
    type = match.arg(type)
      
    # Settings:
    if (is.null(param)) param = evParam(type)$param
    
    # Random Variates:
    ans = revCopula(n = n, param = param, type = type) 
      
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------

    
evCopulaFit =
function(u, v = NULL, type = evList(), ...)
{   # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Fits the paramter of an elliptical copula
    
    # Note:
    #   The upper limit for nu is 100
    
    # FUNCTION:
    
    # Match Arguments:
    type = match.arg(type)
    
    # Settings:
    U <<- u
    V <<- v
    if (is.list(u)) {
        U <<- u[[1]]
        V <<- u[[2]]
    }
    if (is.matrix(u)) {
        U = u[, 1]
        V = u[, 2]
    }

    # Start Values:
    param = evParam(type)$param
    range = evRange(type)
    paramLength = length(param)
    
    # Log-Likelihood Function:
    .fun = function(x, type) {
        -mean( log(devCopula(u = U, v = V, param = x, type = type)) )
    }
        
    if (paramLength == 1) {
        # We have only one parameter to optimize ...
        fit = optimize(f = .fun, lower = range[1], upper = range[2], 
            maximum = FALSE, tol = .Machine$double.eps^0.25, 
            type = type, ...)
    } else {
        # Log-Likelihood Function:
        range = evRange(type)
        fit = nlminb(start = param, objective = .fun, 
            lower = range[1], upper = range[2], type = type, ...)
    }
    
    # Return Value:
    fit
}


################################################################################

