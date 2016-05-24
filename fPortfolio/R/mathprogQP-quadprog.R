
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
#  rquadprogQP                  Rmetrics Interface for QUADPROG QP solvers 
#  quadprogQP                   Convenience wrapper for QUADPROG QP solvers
#  quadprogQPControl            QUADPROG QP control parameter list 
#  rquadprog                    Synonyme name for quadprog::solveLP function        
###############################################################################


rquadprogQP <- 
    function(objective, lower=0, upper=1, linCons, control=list())
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Implements Goldberg-Idnani Algorithm
    
    # Arguments:
    #   objective - list(dvec=NULL, Dmat=NULL)
    #   lower - lower box constraints
    #   upper - upper box constraints
    #   linCons - linear constraints, list with entries:
    #       mat, lower, upper.
    #   control - control list
    
    # FUNCTION:
    
    # Control List:
    ctrl <- quadprogQPControl()
    if (length(control) > 0)
        for (name in names(control)) ctrl[name] = control[name]
    control <- ctrl
    
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
      
    # Linear Constraints:  
    if(missing(linCons)) {
        eqA <- ineqA <- NULL
        eqA.bound <- ineqA.lower <- ineqA.upper <- NULL
    } else {
        mat <- linCons[[1]]
        M <- nrow(mat)
        lower <- linCons[[2]]
        upper <- linCons[[3]]
        if(length(lower) == 1) lower <- rep(lower, M)
        if(length(upper) == 1) upper <- rep(upper, M)
        
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
    optim <- quadprogQP(
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
        solver = "quadprog",
        version = optim$version)
    class(value) = c("solver", "list")
    value
}


###############################################################################


quadprogQP <- 
    function(
        objective=list(dvec=NULL, Dmat=NULL), 
        par.lower=NULL, par.upper=NULL, 
        eqA=NULL, eqA.bound=NULL,  
        ineqA=NULL, ineqA.lower=NULL, ineqA.upper=NULL,
        control=list())
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Implements Goldberg-Idnani Algorithm
    
    # FUNCTION:
    
    # Control List:
    ctrl <- quadprogQPControl()
    if (length(control) > 0)
        for (name in names(control)) ctrl[name] = control[name]
    control <- ctrl
    
    # General Settings:
    dvec <- -objective$dvec
    Dmat <- objective$Dmat
    Names <- colnames(rbind(dvec, Dmat))
    N <- ncol(rbind(dvec, Dmat))
    
    # Box Constraints:
    if (length(par.lower) == 1) par.lower <- rep(par.lower, N)
    if (length(par.upper) == 1) par.upper <- rep(par.upper, N)

    # Constraints Settings:
    Amat <- eqA
    if (!is.null(ineqA)) Amat <- rbind(eqA, ineqA, -ineqA)
    Amat <- rbind(Amat, diag(N), -diag(N))  
    Amat <- t(Amat) 
    bvec <- eqA.bound
    if (!is.null(ineqA.lower)) bvec <- c(bvec, ineqA.lower)
    if (!is.null(ineqA.upper)) bvec <- c(bvec, -ineqA.upper)
    bvec <-  c(bvec, par.lower)
    if (!is.null(par.upper)) bvec <- c(bvec, -par.upper)
    if (is.null(eqA)) meq <- 0 else meq <- nrow(eqA)
    Amat <- Amat[, is.finite(bvec)]
    bvec <- bvec[is.finite(bvec)]
    
    # Optimize:
    elapsed <- Sys.time()
    optim <- quadprog::solve.QP(
        Dmat = Dmat, 
        dvec = dvec, 
        Amat = Amat, 
        bvec = bvec, 
        meq = meq)
    elapsed <- Sys.time() - elapsed
        
    # Note:
    #   DW: if quadprog::solve.QP fails with non-zero status optim$ierr=
    #   =1: stop("constraints are inconsistent, no solution!")
    #   =2: stop("matrix D in quadratic function is not positive definite!")
    #   this is ugly!
    #   Thus:
    Status <- 0
    Message <- "solution found"
    
    # Add:  
    names(optim$solution) <- Names
    
    # Version:
    package <- packageDescription(pkg="quadprog")
    version <- paste(package$Package, package$Version, package$Date)
        
    # Return Value:
    value <- list(
        opt = optim, 
        solution = optim$solution, 
        objective = optim$value,
        status = Status,
        message = Message,
        solver = "quadprog",
        elapsed <- elapsed,
        version = version)
    class(value) = c("solver", "list")
    value
}


###############################################################################


rquadprog <- quadprog::solve.QP


###############################################################################


quadprogQPControl <-
    function(solver="quadprog", trace=FALSE)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns control parameter list
    
    # Arguments:
    #   trace - al logical flag, should the function be traced?
    
    # Details:
    #   Note there are no control paramters supported in 
    #   quadprog::solve.QP
    
    # FUNCTION:
    
    # Control Parameter:
    control <- list(trace = trace)
    
    # Return Value:
    control
}


###############################################################################

