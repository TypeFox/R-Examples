
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


################################################################################
# FUNCTION:                DESCRIPTION:
#  rsolnpNLP                NLP wrapper for solver solnpNLP()
#  solnpNLP                 NLP wrapper for solver solnp()
#  solnpNLPControl          Control parameter list 
# REQUIRES:
#  solnp                    Solver in package Rsolnp
################################################################################


rsolnpNLP <- 
    function(start, objective, lower=0, upper=1, linCons, funCons, 
    control=list())
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Function wrapper for solver solnp()
    
    # FUNCTION:
    
    # Load:
    # require(Rsolnp)
    
    # Update Control List:
    ctrl <- solnpNLPControl()
    if (length(control) > 0)
        for (name in names(control)) ctrl[name] = control[name]
    control <- ctrl
    BIG <- 1e6
    N <- length(start)
    
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
        eqA.bound <- ineqA.lower <- ineqA.upper <- NULL
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
    
    # Nonlinear Constraints:
    if(missing(funCons)) {
        eqFun <- ineqFun <- list()
        eqFun.bound <- ineqFun.lower <- ineqFun.upper <- NULL
    } else {
        fun <- funCons[[1]]
        lower <- funCons[[2]]
        upper <- funCons[[3]]
        eqIndex <- which(lower == upper)
        ineqIndex <- which(lower != upper)
        if (length(eqIndex) == 0) {
            eqFun <- list()
            eqFun.boud <- NULL
        } else {
            eqFun = fun[eqIndex]
            eqFun.bound = lower[eqIndex]
        }
        if (length(ineqIndex) == 0) {
            ineqFun <- list()
            ineqFun.lower <- NULL
            ineqFun.upper <- NULL
        } else {
            ineqFun = fun[ineqIndex]
            ineqFun.lower <- lower[ineqIndex]
            ineqFun.upper <- upper[ineqIndex]
        }
    }
    
    # Optimize Portfolio:
    elapsed <- Sys.time()
    optim <- solnpNLP(
        start = start, 
        objective = objective, 
        par.lower = par.lower, 
        par.upper = par.upper,
        eqA = eqA, 
        eqA.bound = eqA.bound,
        ineqA = ineqA, 
        ineqA.lower = ineqA.lower, 
        ineqA.upper = ineqA.upper,
        eqFun = eqFun, 
        eqFun.bound = eqFun.bound,
        ineqFun = ineqFun, 
        ineqFun.lower = ineqFun.lower, 
        ineqFun.upper = ineqFun.upper,
        control = control) 
    elapsed <- Sys.time() - elapsed
    
    # Version:
    package <- packageDescription(pkg="Rsolnp")
    version <- paste(package$Package, package$Version, package$Date)
          
    # Return Value:
    value <- list(
        opt = optim,
        solution = optim$solution, 
        objective = objective(optim$solution), 
        status = optim$status,
        message = optim $message,
        solver = "solnpNLP",
        elapsed = elapsed,
        version = version)
    class(value) <- c("solver", "list")
    
    # Return Value:
    value
}


################################################################################


solnpNLP <- 
    function(
    start, objective,
    par.lower = NULL, par.upper = NULL,
    eqA = NULL, eqA.bound = NULL,
    ineqA = NULL, ineqA.lower = NULL, ineqA.upper = NULL,
    eqFun = list(), eqFun.bound = NULL,
    ineqFun = list(), ineqFun.lower = NULL, ineqFun.upper = NULL,
    control = list())
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Universal function wrapper for solver solnp(). 
    
    # FUNCTION:
    
    # Load:
    # require(Rsolnp)
    fun <- objective
    
    # Control List:
    ctrl <- solnpNLPControl()
    if (length(control) > 0)
        for (name in names(control)) ctrl[name] = control[name]
    control <- ctrl
    
    # Environment Setting:
    env <- .GlobalEnv
    
    # Box Constraints:
    BIG <- 1e8 # inf does not works here, DW.
    if (is.null(par.lower)) par.lower <- rep(-BIG, length(start))
    if (is.null(par.upper)) par.upper <- rep(+BIG, length(start))
    
    # Linear and Function Equality Constraints:
    if (length(eqA) > 0 || length(eqFun) > 0 ) {
        eqfun <- function(x) {
            ans <- NULL
            if(!is.null(eqA)) ans <- c(ans, as.vector(eqA %*% x))
            if (length(eqFun) > 0)
                for (i in 1:length(eqFun)) ans <- c(ans, eqFun[[i]](x))
            ans }
    } else {
        eqfun <- NULL
    }      
    eqB <- c(eqA.bound, eqFun.bound)
    
    # Linear and Function Inequality Constraints:
    if (length(ineqA) > 0 || length(ineqFun) > 0) {
        ineqfun <- function(x) {
            ans <- NULL
            if(!is.null(ineqA)) ans <- c(ans, as.vector(ineqA %*% x))
            if (length(ineqFun) > 0)
                for (i in 1:length(ineqFun)) ans <- c(ans, ineqFun[[i]](x))
            ans }
    } else {
        ineqfun <- NULL
    }  
    ineqLB <- c(ineqA.lower, ineqFun.lower)
    ineqUB <- c(ineqA.upper, ineqFun.upper)
        
    # Optimize Portfolio:
    elapsed <- Sys.time()
    optim <- Rsolnp::solnp(
        pars = start, 
        fun = fun, 
        eqfun = eqfun, 
        eqB = eqB, 
        ineqfun = ineqfun, 
        ineqLB = ineqLB, 
        ineqUB = ineqUB,  
        LB = par.lower, 
        UB = par.upper, 
        control = control) 
    elapsed <- Sys.time() - elapsed
    names(optim$pars) <- names(start)
    
    # Version:
    package <- packageDescription(pkg="Rsolnp")
    version <- paste(package$Package, package$Version, package$Date)
       
    # Return Value:
    value <- list(
        opt = optim,
        solution = optim$pars, 
        objective = fun(optim$pars), 
        status = optim$convergence,
        message = "not available",
        solver = "solnpNLP",
        elapsed = elapsed,
        version = version)
    class(value) <- c("solver", "list")
    value
}


###############################################################################


solnpNLPControl <-
    function(
    rho = 1, 
    outer.iter = 400, 
    inner.iter = 800, 
    delta = 1.0e-7, 
    tol = 1.0e-8, 
    trace = 0)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns control list
    
    # Arguments:
    #   rho - a numeric value. The penalty parameter 
    #   majit - an integer value. The maximum number of major iterations 
    #   minit - an integer value. The maximum number of minor iterations 
    #   delta - a numeric value. The relative step size in forward 
    #       difference evaluation 
    #   tol - a numeric value. The tolerance on feasibility and optimality 
    
    # Notes:
    #   DW: default trace=1 changed to trace=0
    
    # FUNCTION:
    
    # Control Parameters:
    control <- list(
        rho = rho,
        outer.iter = outer.iter,
        inner.iter = inner.iter,
        delta = delta,
        tol = tol,
        trace = trace)
       
    # Return Value:
    control
}


###############################################################################

