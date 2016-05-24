
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
#  rsolveLP                     General function wrapper for LP solvers
#  .solveLP.MAD.demo            MAD portfolio example
#  .solveLP.GLPK.demo           GLPK help page example
###############################################################################


rsolveLP <- 
  function(objective, lower=0, upper=1, linCons, 
           control=list(solver="glpk", invoke=c("R", "AMPL", "NEOS")))
  {    
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Rmetrics Interface for AMPL/NEOS LP solvers
    
    # Argments:
    #    objective - numeric vector.
    #    lwer, upper - box constraints
    #    linCons - linear constraints: mat, lower and upper
    #    control - control list
    
    # Control:
    solver <- control$solver
    invoke <- control$invoke[1]
    
    # Solve Linear Problem:
    if (invoke == "R") {
      rfooLP <- match.fun ( paste("r", solver, "LP", sep=""))
      ans <- rfooLP(objective, lower, upper, linCons, control)
    }
    if (invoke == "AMPL" ) {
      ans <- ramplLP(objective, lower, upper, linCons, control)
    }
    if (invoke == "NEOS" ) {
      ans <- rneosLP(objective, lower, upper, linCons, control)
    }  
    ans$solver <- paste(invoke, ans$solver) 
    
    # Return Value:
    ans  
  }


###############################################################################


.solveLP.MAD.demo <- 
  function()
  {
    # Solve MAD Portfolio:
    
    # Load Dataset
    dataSet <- data("LPP2005REC", package="timeSeries", envir=environment())
    LPP2005REC <- get(dataSet, envir=environment())
    
    # Load Swiss Pension Fund Data:
    nAssets <- 6
    nScenarios <- 100   
    data <- 100 * LPP2005REC[1:nScenarios, 1:nAssets] 
    Mean <- colMeans(data)
    Data <- as.matrix(data)
    targetReturn <- mean(data)
    
    # Objective Function:
    vec <- c(weights=rep(0, nAssets), scenarios=rep(1/nScenarios, nScenarios))
    
    # Set up Constraints Matrix:
    mat <- rbind(
      MAD.LE = cbind(Data, -diag(nScenarios)),
      MAD.GE = cbind(Data, +diag(nScenarios)),
      RETURN = t(c(Mean, rep(0, nScenarios))),
      BUDGET = t(c(rep(1, nAssets), rep(0, nScenarios))),
      X = cbind(matrix(rep(0, nAssets*nScenarios),ncol=nAssets), diag(nScenarios)),
      WEIGHTS = cbind(diag(nAssets), matrix(rep(0, nScenarios*nAssets), nrow=nAssets)))
    
    # Set up Right Hand Side of Constraints Equations:
    rhs <- c(
      MAD.LE = rep(0, nScenarios),
      MAD.GE = rep(0, nScenarios),
      RETURN = targetReturn,
      BUDGET = 1,
      X = rep(0, nScenarios),
      WEIGHTS = rep(0, nAssets))
    
    # Set up Vector of Directions:
    dir <- c(
      MAD.LE = rep("<=", nScenarios),
      MAS.GE = rep(">=", nScenarios),
      RETURN = "==",
      BUDGET = "==",
      X = rep(">=", nScenarios),
      WEIGHTS = rep(">=", nAssets))
    
    # Conversions:
    RHS <- rep(Inf, times=length(dir))
    RHS[dir == "<="] <- rhs[dir == "<="]
    RHS[dir == "=="] <- rhs[dir == "=="] 
    LHS <- rep(-Inf, times=length(dir))
    LHS[dir == ">="] <- rhs[dir == ">="]
    LHS[dir == "=="] <- rhs[dir == "=="]
    
    # Arguments: 
    objective <- vec
    lower <- 0
    upper <- 1
    linCons <- list(mat, LHS, RHS)
    control <- list()
    
    # Contributed R Solver - Original Function Calls:
    Rglpk::Rglpk_solve_LP(vec, mat, dir, rhs)
    
    # Contributed R Solver - Interfaced  
    rglpkLP(objective, lower, upper, linCons)
    
    # AMPL:
    ramplLP(objective, lower, upper, linCons)
    
    # All AMPL:
    for (solver in c(
      "cplex", "donlp2", "gurobi", "loqo", "lpsolve", "minos", "snopt",
      "ipopt", "bonmin", "couenne")) {
      ans <- ramplLP(objective=vec, lower, upper, linCons, control=list(solver=solver))
      print(ans)
    }
    
    # NEOS:
    # require(rneos)
    # lp: Linear Programming Solver:
    for (solver in c("gurobi", "mosek", "ooqp"))   
      print(rneosLP(objective=vec, lower, upper, linCons, 
                    control=list(solver=solver, category="lp"))) 
    # nco: Using Nonlinear Constrained Optimization Solver:    
    for (solver in c(
      "conopt", "filter", "knitro", "lancelot", "loqo", "minos", "mosek",
      "pennon", "snopt"))    
      print(rneosLP(objective=vec, lower, upper, linCons, 
                    control=list(solver=solver, category="nco")))
  }


# -----------------------------------------------------------------------------


.solveLP.GLPK.demo <- 
  function()
  {
    # GLPK Demo from Help Rglpk Page:
    vec <- -c(2, 4, 3)
    mat <- matrix(c(
      3, 2, 1, 
      4, 1, 3, 
      2, 2, 2), 3, 3)
    dir <- c("<=", "<=", "<=")
    rhs <- c(60, 40, 80) 
    
    # For Testing:
    # mat <- rbind(mat, c(0,0,0))
    # dir <- c(dir, "<=")
    # rhs <- c(rhs, 1000)
    
    # Arguments:
    objective <- vec
    lower <- 0
    upper <- Inf
    linCons <- list(mat, lower=-Inf, upper=rhs)
    control <- list()
    
    # Contributed R Solver - Original Function Calls:
    Rglpk::Rglpk_solve_LP(vec, mat, dir, rhs)
    
    # Contributed R Solver - Interfaced  
    rglpkLP(objective, lower, upper, linCons)
    
    # AMPL:
    ramplLP(objective, lower, upper, linCons)
    
    # All AMPL:
    for (solver in c(
      "cplex", "donlp2", "gurobi", "loqo", "lpsolve", "minos", "snopt",
      "ipopt", "bonmin", "couenne")) {
      ans <- ramplLP(objective=vec, lower, upper, linCons, control=list(solver=solver))
      print(ans)
    }
    
    # NEOS:
    # require(rneos)
    # lp: Linear Programming Solver:
    for (solver in c("gurobi", "mosek", "ooqp"))   
      print(rneosLP(objective=vec, lower, upper, linCons, 
                    control=list(solver=solver, category="lp"))) 
    # nco: Using Nonlinear Constrained Optimization Solver:    
    for (solver in c(
      "conopt", "filter", "knitro", "lancelot", "loqo", "minos", "mosek",
      "pennon", "snopt"))    
      print(rneosLP(objective=vec, lower, upper, linCons, 
                    control=list(solver=solver, category="nco")))
  }


###############################################################################

