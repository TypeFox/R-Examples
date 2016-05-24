
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


.exampleData <-
  function()
  {
    
    dataSet <- data("LPP2005REC", package="timeSeries", envir=environment())
    LPP2005REC <- get(dataSet, envir=environment())
    
    nAssets = 6
    lppData = 100 * (as.matrix(LPP2005REC[,-1]))[1:20, 1:nAssets]
    nScenarios = nrow(lppData)
    
    
    targetReturn = mean(lppData)
    Mean = colMeans(lppData)
    Cov = cov(lppData)
    
    start = rep(1/nAssets, nAssets)
    names(start) = colnames(lppData)
    
    madObjective = c( rep(1, nScenarios), rep(0, nAssets)) / nScenarios
    markowitzObjective = list(dvec = rep(0, nAssets), Dmat = Cov)
    sharpeObjective = function(x)  { ( -Mean %*% x / sqrt ( x %*% Cov %*% x )[[1]] ) }
    
    
    markowitzFun = function(x) { ( x %*% Cov %*% x)[[1]] }
    sharpeFun = function(x)  { ( -Mean %*% x / sqrt ( x %*% Cov %*% x )[[1]] ) }
    
    
    lower = 0
    upper = 1
    
    
    # For LP with glpk  
    vec <- c( rep(1, nScenarios), rep(0, nAssets)) / nScenarios
    A <- rbind(
      c( rep(0, nScenarios), targetReturn = Mean),
      c( rep(0, nScenarios), budget = rep(1, nAssets)),
      cbind(-diag(nScenarios), Returns = lppData),
      cbind(+diag(nScenarios), Returns = lppData))
    dirA <- c( 
      rep("==", 2),
      rep("<=", nScenarios),
      rep(">=", nScenarios))      
    rhsA <- c(
      targetReturn = targetReturn, 
      budget = 1,
      lower = rep(0, nScenarios),
      upper = rep(0, nScenarios))
    
    
    # MAD with Box Constraints:
    madbox.linCons <- list(
      mat = rbind(
        c( rep(0, nScenarios), targetReturn = Mean),
        c( rep(0, nScenarios), budget = rep(1, nAssets)),
        cbind(-diag(nScenarios), Returns = lppData),
        cbind(+diag(nScenarios), Returns = lppData)),
      lower = c(
        targetReturn = targetReturn, 
        budget = 1, 
        rep(-Inf, nScenarios), 
        rep(   0, nScenarios)),
      upper = c(
        targetReturn = targetReturn, 
        budget = 1, 
        rep(  0, nScenarios), 
        rep(Inf, nScenarios)))
    
    
    
    # MAD with Box and Group Constraints:
    Bonds    = c(1, 0, 0, 1, 0, 0);  bonds    = 0.3  # >=  30%
    Foreign  = c(0, 0, 0, 1, 1, 1);  foreign  = 0.5  # <=  50%
    Equities = c(0, 0, 1, 0, 1, 1);  equities = 0.6  # <=  60%
    madgroup.linCons = list(
      mat = rbind(
        c( rep(0, nScenarios), targetReturn = Mean),
        c( rep(0, nScenarios), budget = rep(1, nAssets)),
        c( rep(0, nScenarios), bons = Bonds),
        c( rep(0, nScenarios), foreign = Foreign),
        c( rep(0, nScenarios), equities = Equities),
        cbind(negDiag = -diag(nScenarios), Returns = lppData),
        cbind(posDiag = +diag(nScenarios), Returns = lppData)),
      lower = c(
        targetReturn = targetReturn, 
        budget = 1, 
        bonds = 0.3,
        foreign = 0,
        equities = 0,
        rep(-Inf, nScenarios), 
        rep(   0, nScenarios)),
      upper = c(
        targetReturn = targetReturn, 
        budget = 1, 
        bonds = 1,
        foreign = 0.5,
        equities = 0.6,
        rep(   0, nScenarios), 
        rep(Inf, nScenarios)))  
    
    
    # QP Box Constrained Markowitz:
    Budget = rep(1, nAssets)
    investBudget = 1
    qpbox.linCons <- list(
      mat = rbind(Return = Mean, Budget = Budget),
      lower = c(Return=targetReturn, Budget=investBudget), 
      upper = c(Return=targetReturn, Budget=investBudget) )  
    
    
    # QP Box and Group Constrained Markowitz:
    Bonds    = c(1, 0, 0, 1, 0, 0);  bonds    = 0.3  # >=  30%
    Foreign  = c(0, 0, 0, 1, 1, 1);  foreign  = 0.5  # <=  50%
    Equities = c(0, 0, 1, 0, 1, 1);  equities = 0.6  # <=  60%    
    qpgroup.linCons <- list(
      mat = rbind(Return = Mean, Budget = Budget, 
                  Bonds = Bonds, Foreign = Foreign, Equtities = Equities),
      lower = c(Return=targetReturn, Budget=1, 
                Bonds=0.3, Foreign=0.0, Equtities=0.0), 
      upper = c(Return=targetReturn, Budget=1, 
                Bonds=1.0, Foreign=0.5, Equtities=0.6) )  
    
    
    # QP Box Constrained Markowitz:     
    qpbox.funCons <- list(
      fun = list(
        Return = function(x) { Mean %*%  x },
        Budget = function(x) { sum(x) } ),
      lower = c(Return=targetReturn, Budget=1), 
      upper = c(Return=targetReturn, Budget=1) )          
    
    
    
    qpgroup.funCons <- list(
      fun = list(
        Return = function(x) { Mean %*%  x },
        Budget = function(x) { sum(x) },
        Bonds = function(x) { Bonds %*% x },
        Foreign = function(x) { Foreign %*% x },
        Equtities = function(x) { Equities %*% x } ),
      lower = c(Return=targetReturn, Budget=1, Bonds=0.3, Foreign=0.0, Equtities=0.0), 
      upper = c(Return=targetReturn, Budget=1, Bonds=1.0, Foreign=0.5, Equtities=0.6) )        
    
    
    budget.linCons <- list(
      mat = rbind(Budget = Budget),
      lower = investBudget, 
      upper = investBudget)  
    
    budget.funCons <- list(
      fun = list(function(x) sum(x)),
      lower = investBudget, 
      upper = investBudget)  
    
    
    groupbudget.funCons <- list(
      fun = list(
        Budget = function(x) { sum(x) },
        Bonds = function(x) { Bonds %*% x },
        Foreign = function(x) { Foreign %*% x },
        Equtities = function(x) { Equities %*% x } ),
      lower = c(Budget=1, Bonds=0.3, Foreign=0.0, Equtities=0.0), 
      upper = c(Budget=1, Bonds=1.0, Foreign=0.5, Equtities=0.6) )              
    
  }


###############################################################################

