
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
# FUNCTION:                  ARCHIMEDEAN COPULAE PARAMETER FITTING:
#  archmCopulaSim             Simulates bivariate elliptical copula
#  archmCopulaFit             Fits the paramter of an elliptical copula
################################################################################


################################################################################
# FUNCTION:                  ARCHIMEDEAN COPULAE PARAMETER FITTING:
#  archmCopulaSim             Simulates bivariate elliptical copula
#  archmCopulaFit             Fits the paramter of an elliptical copula


archmCopulaSim <-  
    function (n, alpha = NULL, type = archmList()) 
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Simulates bivariate elliptical Copula
    
    # Match Arguments:
    type <- match.arg(type)
    Type <- as.integer(type)
      
    # Settings:
    if (is.null(alpha)) alpha = archmParam(type)$param
    
    # Random Variates:
    ans <- rarchmCopula(n = n, alpha = alpha, type = type) 

    # Control:
    control = list(alpha = alpha[[1]], copula = "archm", type = type)
    attr(ans, "control")<-unlist(control)
        
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------

    
archmCopulaFit <- 
    function(u, v = NULL, type = archmList(), ...)
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Fits the paramter of an elliptical copula
    
    # Note:
    #   The upper limit for nu is 100
    
    # FUNCTION:
    
    # Match Arguments:
    type = match.arg(type)
    Type = as.integer(type)
    
    # Settings:
    U = u
    V = v
    if (is.list(u)) {
        U = u[[1]]
        V = u[[2]]
    }
    if (is.matrix(u)) {
        U = u[, 1]
        V = u[, 2]
    }

    # Estimate Rho from Kendall's tau for all types of Copula:
    alpha = archmParam(type)$param
     
    # Estimate Copula:
    fun = function(x, type, U, V) {
        -mean( log(darchmCopula(u = U, v = V, alpha = x, type = type)) )
    }
    range = archmRange(type)

    fit = nlminb(start = alpha, objective = fun, 
        lower = range[1], upper = range[2],  type = type, U = U, V = V, ...)
      
    # Return Value:
    fit
}


################################################################################

