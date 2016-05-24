
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
# FUNCTION:                    DESCRIPTION:
#  rdonlp2NLP                   Rmetrics Interface for DONLP2 LP sol
#  donlp2NLP                    Convenience wrapper for DONLP2 LP so    
#  donlp2NLPControl             DONLP2 LP control parameter list
#  rdonlp2                      Synonyme name for Rdonlp2::donlp2 function
################################################################################


rdonlp2NLP <- 
    function(
        start, objective, lower=0, upper=1, linCons, funCons, control=list())
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Function wrapper for solver donlp2()
    
    # FUNCTION:
    
    # Update Control List:
    ctrl <- donlp2NLPControl()
    if (length(control) > 0)
        for (name in names(control)) ctrl[name] = control[name]
    control <- ctrl
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
    
    # Linear Constraints:  
    if(missing(linCons)) {
        eqA <- ineqA <- NULL
        eqA.bound <- ineqA.lower <- ineqA.upper <- NULL
    } else {
        mat <- linCons[[1]]
        lower <- linCons[[2]]
        upper <- linCons[[3]]
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
            eqFun <- fun[eqIndex]
            eqFun.bound <- lower[eqIndex]
        }
        if (length(ineqIndex) == 0) {
            ineqFun <- list()
            ineqFun.lower <- NULL
            ineqFun.upper <- NULL
        } else {
            ineqFun <- fun[ineqIndex]
            ineqFun.lower <- lower[ineqIndex]
            ineqFun.upper <- upper[ineqIndex]
        }
    }
    
    # Optimize Portfolio:
    optim <- donlp2NLP(
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
          
    # Return Value:
    value = list(
        opt = optim,
        solution = optim$solution, 
        objective = objective(optim$solution)[[1]], 
        status = optim$status,
        message = optim$message,
        solver = "donlp2NLP")
    class(value) <- c("solver", "list")
    value
}


################################################################################


donlp2NLP <- 
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
    #   NLP wrapper for solver donlp2
    
    # FUNCTION:
    
    # Environment:
    env <- .GlobalEnv
    
    # Update Control List:
    ctrl <- donlp2NLPControl()
    if (length(control) > 0)
        for (name in names(control)) ctrl[name] <- control[name]
    control <- ctrl
    
    # Set Box Constraints:
    if (is.null(par.lower)) par.lower <- rep(-Inf, length(start))
    if (is.null(par.upper)) par.upper <- rep(+Inf, length(start))
    if (length(par.lower) == 1) par.lower <- rep(par.lower, length(start))
    if (length(par.upper) == 1) par.upper <- rep(par.upper, length(start))
    
    # Set Linear Equality and Inequality Constraints:
    A <- rbind(eqA, ineqA)
    lin.lower <- c(eqA.bound, ineqA.lower)
    lin.upper <- c(eqA.bound, ineqA.upper)
    
    # Set Nonlinear Equality and Inequality Constraints:
    if ((length(eqFun) + length(ineqFun)) == 0) {
        nlin <- list()
        nlin.lower <- rep(-Inf, length(nlin))
        nlin.upper <- rep(+Inf, length(nlin))
    } else {
        nlin <- list()
        if (length(eqFun) > 0) nlin = c(nlin, eqFun)
        if (length(ineqFun) > 0) nlin = c(nlin, ineqFun)
        nlin.lower <- c(eqFun.bound, ineqFun.lower)
        nlin.upper <- c(eqFun.bound, ineqFun.upper)
    }
    
    # Optimize Portfolio:
    optim <- rdonlp2(
        par = start, 
        fn = objective, 
        par.upper = par.upper, 
        par.lower = par.lower, 
        A = A, 
        lin.upper = lin.upper, 
        lin.lower = lin.lower, 
        nlin = nlin, 
        nlin.upper = nlin.upper, 
        nlin.lower = nlin.lower, 
        control = control, 
        control.fun = function(lst) {return(TRUE)}, 
        env = .GlobalEnv, 
        name = NULL)
    names(optim$par) <- names(start) 
    
    # Extract Weights:
    weights <- .checkWeights(optim$par)
    attr(weights, "invest") <- sum(weights)
    
    # Check Messages and Get Status:
    #   ... unfortunately donlp2 has no status vaqriable, 
    #       so we have to analyze the messages
    Status <- 1
    # Message <- "1234567890123456789012345"
    message11 <- "KT-conditions satisfied, " # no further correction computed"
    message12 <- "computed correction small" # , regular case"
    message13 <- "stepsizeselection: x almo" # st feasible, dir. deriv. very small"
    if (substr(optim$message, 1, 25) == message11) Status <- 0 
    if (substr(optim$message, 1, 25) == message12) Status <- 0 
    if (substr(optim$message, 1, 25) == message13) Status <- 0 
       
    # Return Value:
    value <- list(
        opt = optim,
        solution = optim$par, 
        objective = objective(optim$par)[[1]], 
        status = Status,
        message = optim$message,
        solver = "donlp2NLP")
    class(value) <- c("solver", "list")
    value
}


################################################################################


rdonlp2 <- function(...) {
  Rdonlp2::donlp2(...)
}


################################################################################


donlp2NLPControl <-
function (
    iterma = 4000, nstep = 20, fnscale = 1, report = FALSE, 
    rep.freq = 1, tau0 = 1, tau = 0.1, del0 = 1, epsx = 1e-05, 
    delmin = 0.1 * del0, epsdif = 1e-08, nreset.multiplier = 1, 
    difftype = 3, epsfcn = 1e-16, taubnd = 1, hessian = FALSE, 
    te0 = TRUE, te1 = FALSE, te2 = FALSE, te3 = FALSE, silent = TRUE, 
    intakt = TRUE) 
{
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION:
    
    # Control Parameters:
    control <- list(
        iterma = iterma, nstep = nstep, fnscale = fnscale, 
        report = report, rep.freq = rep.freq, tau0 = tau0, tau = tau, 
        del0 = del0, epsx = 1e-05, delmin = delmin, epsdif = epsdif, 
        nreset.multiplier = nreset.multiplier, difftype = difftype, 
        epsfcn = epsfcn, taubnd = taubnd, hessian = hessian, 
        te0 = te0, te1 = te1, te2 = te2, te3 = te3, silent = silent, 
        intakt = intakt)
        
    # Return Value:
    control
}


################################################################################

