
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR Description. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA 02111-1307 USA


################################################################################
# FUNCTION:                DESCRIPTION:
#  solveRdonlp2             Portfolio interface to solver Rdonlp2
#  .rdonlp2Arguments        Returns arguments for solver
################################################################################


solveRdonlp2 <-
  function(data, spec, constraints)
  {
    # Description:
    #   Portfolio interface to solver Rdonlp2
    
    # Arguments;
    #   data - an object of class timeSeries
    #   spec - an object of class fPFOLIOSPEC
    #   constraints - an object of class character
    
    # Example:
    #   
    
    # FUNCTION:   
    
    # Settings:
    Data <- portfolioData(data, spec)
    nAssets <- getNAssets(Data)
    mu <- getMu(Data)
    Sigma <- getSigma(Data)
    
    # Compose Arguments for Solver:
    args <- .rdonlp2Arguments(data, spec, constraints)
    
    # Solve Multiassets Portfolio:
    ans <- Rdonlp2::donlp2(
      par = args$par,
      fn = args$fn,
      par.lower = args$par.lower, 
      par.upper = args$par.upper, 
      A = args$A, 
      lin.lower = args$lin.lower, 
      lin.upper = args$lin.upper,
      nlin = args$nlin, 
      nlin.lower = args$nlin.lower, 
      nlin.upper = args$nlin.upper,
      control = Rdonlp2::donlp2Control(),
      control.fun = function(lst) {return(TRUE)}, 
      env = .GlobalEnv, 
      name = NULL)      
    returnFun <- match.fun(getObjective(spec)[2])
    ans$targetReturn <- returnFun(ans$par)
    riskFun <- match.fun(getObjective(spec)[3])
    ans$targetRisk <- riskFun(ans$par)
    ans$ans <- ans
    
    ans$solver <- "solveRdonlp2"
    ans$objective <- args$fn(ans$par)
    names(ans$par) <- names(start)
    ans$weights <- .checkWeights(ans$par)
    ans$solution <- ans$weights
    attr(ans$weights, "invest") <- sum(ans$weights)
    ans$status <- 1
    message11 <- "KT-conditions satisfied, "
    message12 <- "computed correction small"
    message13 <- "stepsizeselection: x almo"
    if (substr(ans$message, 1, 25) == message11) 
      ans$status <- 0
    if (substr(ans$message, 1, 25) == message12) 
      ans$status <- 0
    if (substr(ans$message, 1, 25) == message13) 
      ans$status <- 0
    
    # Return Value:
    class(ans) <- c("solver", "list")
    ans
  }


# ------------------------------------------------------------------------------


.rdonlp2Arguments <-
  function(data, spec, constraints)
  {
    # Description:
    #   Create Arguments for Rdonlp2
    
    # Details:
    #       min:                    fn(x)
    #       subject to: 
    #                     par.lower <= x <= par.upper
    #                  lin.lower <= A %*% x <= lin.upper
    #                 nlin.lower <= nlin(x) <= nlin.upper  
    
    # FUNCTION:  
    
    DEBUG = FALSE
    
    # Settings:
    Data <- portfolioData(data)
    nAssets <- getNAssets(Data)  
    mu <- getMu(Data)
    Sigma <- getSigma(Data)
    fn <- match.fun(getObjective(spec)[1])
    
    # Box Constrains:
    par.lower <- minWConstraints(data, spec, constraints)
    par.upper <- maxWConstraints(data, spec, constraints)
    if(DEBUG) print(rbind(par.lower, par.upper))
    
    # Linear / Group Constraints:
    # ... targetReturn may be not defined,then set it to NA
    if (is.null(getTargetReturn(spec))) setTargetReturn(spec) <- NA
    # ... has in the first line the return constraint, if NA then ignore  it
    eqsumW <- eqsumWConstraints(data, spec, constraints)
    if (is.na(eqsumW[1, 1])) eqsumW = eqsumW[-1, , drop= FALSE]
    Aeqsum <- eqsumW[, -1]
    aeqsum <- eqsumW[, 1]
    minsumW <- minsumWConstraints(data, spec, constraints)
    if (is.null(minsumW)) {
      Aminsum <- aminsum <- NULL
    } else {
      Aminsum <- minsumW[, -1]
      aminsum <- minsumW[, 1]
    }      
    maxsumW <- maxsumWConstraints(data, spec, constraints)
    if (is.null(maxsumW)) {
      Amaxsum <- amaxsum <- NULL
    } else {
      Amaxsum <- maxsumW[, -1]
      amaxsum <- maxsumW[, 1]
    }      
    A <- rbind(Aeqsum, Aminsum, Amaxsum)
    lin.lower <- c(aeqsum, aminsum, rep(-Inf, length(amaxsum)))
    lin.upper <- c(aeqsum, rep(Inf, length(aminsum)), amaxsum)
    if(DEBUG) print(cbind(lin.lower, A, lin.upper))
    
    # Nonlinear Constraints - Here Covariance Risk Budgets:
    nlin <- list()
    nlin.lower <- NULL
    nlin.upper <- NULL
    # Check Constraints Strings for Risk Budgets:
    # Example: constraints = c("minB[2:3]=0.1", "maxB[3:5]=0.9")
    validStrings <- c("minB", "maxB")
    usedStrings <- unique(sort(sub("\\[.*", "", constraints)))
    checkStrings <- sum(usedStrings %in% validStrings)
    includeRiskBudgeting <- as.logical(checkStrings)
    if (DEBUG) print(includeRiskBudgeting)
    if (includeRiskBudgeting) {
      # Compose Non-Linear (Cov Risk Budget) Constraints Functions:
      nlcon <- function(x) {
        B1 <- as.vector(x %*% Sigma %*% x)
        B2 <- as.vector(x * Sigma %*% x)
        B <- B2/B1
        B
      }
      if(DEBUG) print(nlcon)
      # Compose non-linear functions now for each asset ...
      for (I in 1:nAssets)
        eval( parse(text = paste(
          "nlcon", I, " = function(x) { nlcon(x)[", I, "] }", sep = "")) )
      nlinFunctions <- paste("nlcon", 1:nAssets, sep = "", collapse = ",")
      nlinFunctions <- paste("list(", nlinFunctions, ")")
      nlin <- eval( parse(text = nlinFunctions) )
      if(DEBUG) print(nlin)
      # ... and finally Compose Constraints Vectors:
      nlin.lower <- minBConstraints(data, spec, constraints)
      nlin.upper <- maxBConstraints(data, spec, constraints)
      if(DEBUG) print(rbind(nlin.lower, nlin.upper))
    }
    
    # General non-lin Portfolio Constraints:
    # ... todo: currently overwrites previous selection
    nlin <- listFConstraints(data, spec, constraints)
    if(DEBUG) print(nlin)
    nlin.lower <- minFConstraints(data, spec, constraints)
    nlin.upper <- maxFConstraints(data, spec, constraints)
    if(DEBUG) print(cbind(nlin.lower, nlin.upper))
    
    # Return Value:
    list(
      par = rep(1/nAssets, nAssets), fn = fn,
      par.lower = par.lower, par.upper = par.upper, 
      A = A, lin.lower = lin.lower, lin.upper = lin.upper,
      nlin = nlin, nlin.lower = nlin.lower, nlin.upper = nlin.upper)
  }


################################################################################

