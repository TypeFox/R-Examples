
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
# FUNCTION:               DESCRIPTION:
#  solveRampl.MV           AMPL solver for a MV Long Only Portfolio
# FUNCTION:               DESCRIPTION:
#  solveRampl.CVAR         AMPL solver for a CVAR Long Only Portfolio
################################################################################


solveRampl.MV <-
    function (data, spec, constraints="LongOnly") 
{
    # A function written by Diethelm Wuertz
    
    # Description:
    #   A Mean-Variance (MV) Long Only portfolio solver
    
    # Note:
    #   No other constraints are allowed.
    
    # Example:
    #    data <- 100 * LPP2005.RET[, 1:6]
    #    spec <- portfolioSpec()
    #    setTargetReturn(spec) <- mean(data)
    #    setSolver(spec) <- "solveRampl.MV"
    #    solveRampl.MV(data, spec)
    
    # FUNCTION:
    
    # Check Consistency:
    stopifnot(getType(spec) == "MV")
    stopifnot(constraints == "LongOnly")
    
    # Force AMPL Settings:
    spec@ampl$ampl <- TRUE
    solver <- spec@ampl$solver
    project <- spec@ampl$project
    trace <- spec@ampl$trace
    
    # Get Portfolio Settings:
    Data <- portfolioData(data, spec)
    Sigma <- getSigma(Data)
    n <- nAssets <- getNAssets(Data)
    targetReturn <- getTargetReturn(spec)    
    
    # Create AMPL Model File:
    amplModelOpen(project)
    model <- c(
        "# Quadratic Programming",
        "# Long Only Markowitz",
        "param nAssets;",                        
        "param mu{1..nAssets};",                   
        "param Sigma{1..nAssets, 1..nAssets};",      
        "param targetReturn;",                                                    
        "var x{1..nAssets} >= 0;",            
        "minimize Risk: sum {i in 1..nAssets} sum{j in 1..nAssets} x[i]*Sigma[i,j]*x[j];",
        "subject to Return: sum{i in 1..nAssets} mu[i]*x[i] = targetReturn;",
        "subject to Budget: sum{i in 1..nAssets} x[i] = 1;",
        NULL )
    amplModelAdd(model, project)
    
    # Create AMPL Data File:
    amplDataOpen(project)
    amplDataAddValue(data="nAssets", value=n, project)
    amplDataAddVector(data="mu", vector=getMu(Data), project)
    amplDataAddMatrix(data="Sigma", matrix=getSigma(Data), project)
    amplDataAddValue(data="targetReturn", value=getTargetReturn(spec), project)
    
    # Create AMPL Run File:
    amplRunOpen(project)
    run <- c(
        paste("reset ;"),
        paste("option solver", solver, ";", sp = ""),
        paste("model ", project, ".mod ;", sep = ""),
        paste("data ", project, ".dat ;", sep = ""),
        paste("solve ;"),
        paste("display x > ", project, ".txt ;", sep = ""),
        paste("display solve_result_num > ", project, ".txt ;", sep = ""),
        paste("display solve_result > ", project, ".txt ;", sep = ""),
        paste("display solve_message > ", project, ".txt ;", sep = ""),
        paste("exit ;"),
        NULL )
    amplRunAdd(run, project)
   
    # Exec AMPL:
    command <- paste("ampl -t -vs", paste(project, "run", sep="."))
    solve <- system(command, intern=TRUE)                          
    
    # Read AMPL Output File: 
    file <- paste(project, "txt", sep = ".")
    out <- scan(file, what = character(0), sep="\n", quiet=TRUE)
  
    # Get Weights:
    Index <- (grep(";", out) - 1)[1]
    splits <- strsplit(paste(out[2:Index], collapse=" "), " ")[[1]]
    solution <- as.numeric(splits[splits != ""])[seq(2, 2*n, by=2)]
    Index <- as.numeric(splits[splits != ""])[seq(1, 2*n, by=2)]
    solution[Index] <- solution
    
    # Get Status:
    status <- strsplit(out[grep("solve_result", out)], split=" ")
    statusCode <- status[[1]][3]
    statusMessage <- status[[2]][3]
    
    # Version:
    version <- system(paste(solver, "-v"), intern=TRUE) 
    
    # Compose Result into a List:
    ans <- list(
        type = getType(spec), 
        solver = getSolver(spec), 
        optim = solve, 
        weights = solution, 
        targetReturn = targetReturn, 
        targetRisk = sqrt(solution %*% Sigma %*% solution)[[1, 1]], 
        objective = sqrt(solution %*% Sigma %*% solution)[[1, 1]], 
        status = statusCode, 
        message = statusMessage)
    class(ans) <- c("solver", "list")
      
    # Return Value: 
    ans
}


################################################################################


solveRampl.CVAR <-
    function(data, spec, constraints="LongOnly")
{    
    # A function written by Diethelm Wuertz
    
    # FUNCTION:
    
    # Portfolio Model
    #   max     CVaR
    #   s.t.    desired target Return
    #   s.t.    full Investement
    #   s.t.    long only positions

    # Check Consistency:
    stopifnot(getType(spec) == "CVAR")
    stopifnot(getSolver(spec) == "solveAMPLipoptCVAR")
    stopifnot(constraints == "LongOnly")
    
    # Force AMPL Settings:
    spec@ampl$ampl <- TRUE
    project <- spec@ampl$project
    trace <- spec@ampl$trace
    
    # Settings:
    Data <- portfolioData(data, spec)
    data <- getSeries(Data)
    n <- nAssets <- getNAssets(Data)
    nScenarios <- nrow(getSeries(Data))
    Mean <- getMean(Data)
    targetReturn <- getTargetReturn(spec)
    alpha <- getAlpha(spec)
    Type <- getType(spec)
    project <- "ampl"
    solver <- "ipopt"

    # Model File:
    amplModelOpen(project)
    model <- c(
        "param alpha ;",
        "param nAssets ;",
        "param nScenarios ;",
        "param Mean{1..nAssets} ;",
        "param targetReturn ;",
        "param Data{1..nScenarios,1..nAssets} ;",
        "var weights{1..nAssets} ;",
        "var VaR ;",
        "var z{1..nScenarios};",
        "maximize Risk: VaR - ( sum{i in 1..nScenarios} z[i] ) / ( alpha  * nScenarios ) ;",
        "subject to Return: sum{i in 1..nAssets} weights[i] * Mean[i] = targetReturn ;" ,
        "subject to Weights{i in 1..nAssets}:  weights[i] >= 0 ;", 
        "subject to Budget: sum{i in 1..nAssets} weights[i] = 1 ;", 
        "subject to Scenarios{k in 1..nScenarios}: -VaR + z[k] + sum{i in 1..nAssets} Data[k,i]*weights[i] >= 0 ;" ,
        "subject to Z{i in 1..nScenarios}:  z[i] >= 0 ;", 
        NULL )
    amplModelAdd(model, project)

    # Data File:
    amplDataOpen(project)
    amplDataAddValue(data="alpha", value=alpha, project)
    amplDataAddValue(data="nAssets", value=nAssets, project)
    amplDataAddValue(data="nScenarios", value=nScenarios, project)
    amplDataAddValue(data="targetReturn", value=targetReturn, project)
    amplDataAddVector(data="Mean", vector=Mean, project)
    amplDataAddMatrix(data="Data", matrix=data, project)

    # Run File:
    amplRunOpen(project)
    run <- c(
        paste("option solver cplex ;"),
        paste("model ", project, ".mod ;", sep = ""),
        paste("data ", project, ".dat ;", sep = ""),
        "solve ;",
        paste("display weights > ", project, ".txt ;", sep = ""),
        paste("display VaR > ", project, ".txt ;", sep = ""),
        paste("exit ;"),
        NULL)
    amplRunAdd(run, project)

    # Exec AMPL:
    command <- paste("ampl -t -vs", paste(project, "run", sep="."))
    solve <- system(command, intern=TRUE)                          
    
    # Read AMPL Output File: 
    file <- paste(project, "txt", sep = ".")
    out <- scan(file, what = character(0), sep="\n", quiet=TRUE)
  
    # Get Weights:
    Index <- (grep(";", out) - 1)[1]
    splits <- strsplit(paste(out[2:Index], collapse=" "), " ")[[1]]
    solution <- as.numeric(splits[splits != ""])[seq(2, 2*n, by=2)]
    Index <- as.numeric(splits[splits != ""])[seq(1, 2*n, by=2)]
    solution[Index] <- solution
    
    # Get Status:
    status <- strsplit(out[grep("solve_result", out)], split=" ")
    statusCode <- status[[1]][3]
    statusMessage <- status[[2]][3]
    
    # Version:
    version <- system(paste(solver, "-v"), intern=TRUE) 

    # Result:
    ans <- list(
        type = "CVAR", 
        solver = "solveRampl.CVAR", 
        optim = optim, 
        weights = weights,
        targetReturn = targetReturn, 
        targetRisk = -optim$optimum,
        objective = -optim$optimum, 
        status = optim$status[[1]],
        message = "")
    class(ans) <- c("solver", "list")
      
    # Return Value: 
    ans 
}


# #############################################################################

