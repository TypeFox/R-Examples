
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
#  rsolveQP                     General Interface for QP solvers
#  .solveQP.MV.demo
###############################################################################


rsolveQP <- 
  function(objective, lower=0, upper=1, linCons, 
           control=list(solver="quadprog", invoke=c("R", "AMPL", "NEOS")))
  {    
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Implements general function wrapper for QP solvers
    
    # Arguments:
    #   objective - list(dvec=NULL, Dmat=NULL)
    #   lower - lower box constraints
    #   upper - upper box constraints
    #   linCons - linear constraints, list with entries:
    #       mat, lower, upper.
    #   control - control list
    
    # FUNCTION:
    
    # Control:
    solver <- control$solver
    invoke <- control$invoke[1]
    
    # Solve Linear Problem:
    if (invoke == "R") {
      rfooLP <- match.fun ( paste("r", solver, "QP", sep=""))
      ans <- rfooLP(objective, lower, upper, linCons, control)
    }
    if (invoke == "AMPL" ) {
      ans <- ramplQP(objective, lower, upper, linCons, control)
    }
    if (invoke == "NEOS" ) {
      ans <- rneosQP(objective, lower, upper, linCons, control)
    }  
    ans$solver <- paste(invoke, ans$solver) 
    
    # Return Value:
    ans  
  }


###############################################################################


.solveQP.MV.demo <- 
  function()
  {
    # Solve Mean-Variance Portfolio:
    
    # Load Dataset
    dataSet <- data("LPP2005REC", package="timeSeries", envir=environment())
    LPP2005REC <- get(dataSet, envir=environment())
    
    # Load Swiss Pension Fund Data:
    nAssets <- 6
    data <- 100 * LPP2005REC[, 1:nAssets]     
    
    # Arguments: 
    objective <- list(dvec=rep(0, nAssets), Dmat=cov(data))
    lower <- 0
    upper <- 1
    mat <- rbind(
      budget = rep(1, times=nAssets), 
      returns = colMeans(data))
    matLower <- c(
      budget = 1, 
      return = mean(data))
    matUpper <- matLower
    linCons <- list(mat, matLower, matUpper)
    control <- list()
    
    # R Contributed Solvers:
    rquadprogQP(objective, lower, upper, linCons)
    ripopQP(objective, lower, upper, linCons)
    
    # Default - AMPL Interface:
    ampl <- ramplQP(objective, lower, upper, linCons)
    ampl
    
    # All AMPL:
    for (solver in c(
      "cplex", "donlp2", "loqo", "lpsolve", "minos", "snopt",
      "ipopt", "bonmin", "couenne")) {
      ans <- ramplQP(objective, lower, upper, linCons, 
                     control=list(solver=solver))
      print(ans)
    }
    
    # NEOS:
    # require(rneos)
    neos <- rneosQP(objective, lower, upper, linCons, 
                    control=list(solver="ipopt", category="nco"))
    neos
    
    # nco: Using Nonlinear Constrained Optimization Solver:    
    for (solver in c(
      "conopt", "filter", "knitro", "lancelot", "loqo", "minos", "mosek",
      "pennon", "snopt")) 
    {   
      ans <- rneosQP(objective, lower, upper, linCons, 
                     control=list(solver=solver, category="nco"))
      print(ans)
    }
    
    # KRESTREL:
    kestrel <- rkestrelQP(objective, lower, upper, linCons, 
                          control=list(solver="loqo"))
    kestrel
  }


###############################################################################

