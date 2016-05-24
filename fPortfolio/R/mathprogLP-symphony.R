
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


###############################################################################
# FUNCTION:                    DESCRIPTION:
#  rsymphonyLP                  Rmetrics Interface for SYMPHONY LP solvers
#  symphonyLP                   Convenience wrapper for SYMPHONY LP solvers
#  symphonyLPControl            SYMPHONY LP control parameter list
###############################################################################


rsymphonyLP <- 
    function(objective, lower=0, upper=1, linCons, control=list())
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Function wrapper for symphony solver 
    
    # Argments:
    #    objective - numeric vector.
    #    lwer, upper - box constraints
    #    linCons - linear constraints: mat, lower and upper
    #    control - control list
    
    # FUNCTION:
    
    # Update Control List:
    ctrl <- symphonyLPControl()
    if (length(control) > 0)
        for (name in names(control)) ctrl[name] = control[name]
    control <- ctrl
    
    # Box Constraints:
    N = length(objective)
    if(length(lower) == 1) {
        par.lower <- rep(lower, N)
    } else {
        par.lower <- lower
    }
    if(length(upper) == 1) {
        par.upper <- rep(upper, N)
    } else {
        par.upper <- upper
    }
    bounds <- list(
        lower = list(ind = 1:N, val = par.lower),
        upper = list(ind = 1:N, val = par.upper))
        
    # Linear Constraints:
    mat <- linCons[[1]]
    M <- nrow(mat)
    lower <- as.vector(linCons[[2]])
    upper <- as.vector(linCons[[3]])
    
    if(length(lower) == 1) {
        lower <- rep(lower, M)
    } else {
        lower <- lower
    }
    if(length(upper) == 1) {
        upper <- rep(upper, M)
    } else {
        upper <- upper
    }
    
    eqIndex <- which(lower == upper)
    ineqIndex <- which(lower != upper)
    
    eqA <- mat[eqIndex, ]
    ineqA <- mat[ineqIndex, ]
    mat <- rbind(eqA, ineqA, ineqA)
    dir <- c(rep("==", length(eqIndex)), rep("<=", length(ineqIndex)), rep(">=", length(ineqIndex)))
    rhs <- c(upper[eqIndex], upper[ineqIndex], lower[ineqIndex])
    mat <- mat[is.finite(rhs), ]
    dir <- dir[is.finite(rhs)]
    rhs <- rhs[is.finite(rhs)]

    # Optimize Portfolio:
    optim <- symphonyLP(
        obj = objective, 
        mat = mat,
        dir = dir,
        rhs = rhs,
        bounds = bounds,
        types = NULL,
        max = FALSE) 
          
    # Version:
    version <- paste(packageDescription("Rsymphony")[1:3], collapse=" ")
    
    # Return Value:
    value = list(
        opt = optim,
        solution = optim$solution, 
        objective = optim$objval, 
        status = optim$status[[1]],
        message = names(optim$status),
        solver = paste("R", control$solver),
        version = version)
    class(value) <- c("solver", "list")
    value
}


# -----------------------------------------------------------------------------


symphonyLP <- Rsymphony::Rsymphony_solve_LP


# -----------------------------------------------------------------------------


symphonyLPControl <- 
    function(solver="symphony", project="r", trace=FALSE) 
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns control parameter list
    
    # FUNCTION:
    
    # Return Value:
    list(solver=solver, trace=trace)
}


###############################################################################

    