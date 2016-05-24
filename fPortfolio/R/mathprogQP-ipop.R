
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
#  ripopQP                      Rmetrics Interface for LOQO QP solver 
#  ipopQP                       Convenience wrapper for LOQO QP solver
#  ipopQPControl                LOQO QP control parameter list   
#  ripop                        Synonyme name for kernlab::ipop function           
###############################################################################

    
ripopQP <- 
    function(objective, lower=0, upper=1, linCons, control=list())
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Implements Vandenberg's LOQO Algorithm
    
    # Arguments:
    #   objective - list(dvec=NULL, Dmat=NULL)
    #   lower - lower box constraints
    #   upper - upper box constraints
    #   linCons - linear constraints, list with entries:
    #       mat, lower, upper.
    #   control - control list
    
    # FUNCTION:
    
    # Control List:
    ctrl <- ipopQPControl()
    if (length(control) > 0)
        for (name in names(control)) ctrl[name] <- control[name]
    control <- ctrl
    BIG <- control$inf
    
    # General Settings:
    dvec <- objective$dvec
    Dmat <- objective$Dmat
    Names <- colnames(rbind(dvec, Dmat))
    N <- ncol(rbind(dvec, Dmat))
    
    # Box Constraints:
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
    par.lower[is.infinite(par.lower)] <- BIG*sign(par.lower[is.infinite(par.lower)])
    par.upper[is.infinite(par.upper)] <- BIG*sign(par.upper[is.infinite(par.upper)])
      
    # Linear Constraints:  
    if(missing(linCons)) {
        eqA <- ineqA <- NULL
        eqA.bound < ineqA.lower <- ineqA.upper <- NULL
    } else {
        mat <- linCons[[1]]
        M <- nrow(mat)
        lower <- linCons[[2]]
        upper <- linCons[[3]]
        if(length(lower) == 1) lower <- rep(lower, M)
        if(length(upper) == 1) upper <- rep(upper, M)
        lower[is.infinite(lower)] <- BIG*sign(lower[is.infinite(lower)])
        upper[is.infinite(upper)] <- BIG*sign(upper[is.infinite(upper)])
        
        eqIndex <- which(lower == upper)
        ineqIndex <- which(lower != upper)
        
        if (length(eqIndex) == 0) {
            eqA <- NULL
            eqA.bound <- NULL
        } else {
            eqA <- mat[eqIndex, ]
            eqA.bound <- lower[eqIndex]
        }
        if (length(ineqIndex) == 0) {
            ineqA <- NULL
            ineqA.lower <- NULL
            ineqA.upper <- NULL
        } else {
            ineqA <- mat[ineqIndex, ]
            ineqA.lower <- lower[ineqIndex]
            ineqA.upper <- upper[ineqIndex]
        }
    }   
    
    # Optimize Portfolio:
    optim <- ipopQP(
        objective, 
        par.lower, par.upper, 
        eqA, eqA.bound,  
        ineqA, ineqA.lower, ineqA.upper,
        control)
    
    # Return Value:
    value <- list(
        opt = optim, 
        solution = optim$solution, 
        objective = optim$objective,
        status = optim$status,
        message = optim$message,
        solver = "ipop",
        version = optim$version)
    class(value) = c("solver", "list")
    value
}


###############################################################################


ipopQP <- 
    function(
        objective=list(dvec=NULL, Dmat = NULL), 
        par.lower=NULL, par.upper=NULL, 
        eqA=NULL, eqA.bound=NULL,
        ineqA=NULL, ineqA.lower=NULL, ineqA.upper=NULL,
        control=list())
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Implements Vandenberg's LOQO Algorithm
    
    # FUNCTION:
    
    # Control List:
    ctrl <- ipopQPControl()
    if (length(control) > 0)
        for (name in names(control)) ctrl[name] <- control[name]
    control <- ctrl
    
    # General Settings:
    dvec <- objective$dvec
    Dmat <- objective$Dmat 
    Names <- colnames(rbind(dvec, Dmat))
    N <- ncol(rbind(dvec, Dmat))
        
    # Solve:
    optim <- ripop(
        c = objective$dvec, 
        H = objective$Dmat, 
        A = rbind(eqA, ineqA), 
        b = c(eqA.bound, ineqA.lower), 
        l = par.lower, 
        u = par.upper, 
        r = c(eqA.bound, ineqA.upper) - c(eqA.bound, ineqA.lower), 
        sigf = control$sigf, 
        maxiter = control$maxiter, 
        margin = control$margin, 
        bound = control$bound, 
        verb = control$verb)
    
    # Add:  
    if (optim@how == "converged") Status <- 0 else Status <- 1
    par <- optim@primal
    names(par) <- Names
    Message <- optim@how 
    
    # Version:
    package <- packageDescription(pkg="kernlab")
    version <- paste(package$Package, package$Version, package$Date)
        
    # Return Value:
    ans <- list(
        opt = optim, 
        solution = par, 
        objective = c( dvec %*% par + 0.5 * par %*% Dmat %*% par )[[1]],
        status = Status,
        message = Message,
        solver = "ipop",
        version = version)
    class(ans) = c("solver", "list")
    ans
}


###############################################################################


ripop <- kernlab::ipop


###############################################################################


ipopQPControl <- 
    function(
        sigf=12, maxiter=400, margin=0.05, bound=10, verb=0, 
        inf=1e12, solver="ipop", trace=FALSE)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns control parameter list
    
    # FUNCTION:
    
    # Control Parameter:
    control <- list(
        sigf = 12, 
        maxiter = 400, 
        margin = 0.05, 
        bound = 10, 
        verb = 0,
        inf = inf,
        solver = solver,
        trace = trace)
        
    # Return Value:
    control
}   
    

###############################################################################

