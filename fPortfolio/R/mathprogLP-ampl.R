
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
#  ramplLP                      Rmetrics Interface for AMPL LP solvers
#  amplLP                       Convenience wrapper for AMPL LP solvers
#  amplLPControl                AMPL LP control parameter list
###############################################################################


ramplLP <- 
    function(objective, lower=0, upper=1, linCons, control=list())
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Implements AMPL LP Interface
    
    # Arguments:
    #   objective - vec
    
    # FUNCTION:
    
    # Control List:
    ctrl <- amplLPControl()
    if (length(control) > 0) 
        for (name in names(control)) ctrl[name] = control[name]
    control <- ctrl
    
    # Controls:
    solver <- control$solver
    project <- control$project
    inf <- control$inf
    trace <- control$trace
    
    # Objective:
    vec <- objective
    
    # Box Constraints:
    replicate <- function(x, n) if(length(x) == 1) rep(x, n) else x
    n <- length(vec)
    x_L <- replicate(lower, n)
    x_U <- replicate(upper, n)
    x_L[is.infinite(x_L)] <- inf*sign(x_L[is.infinite(x_L)])
    x_U[is.infinite(x_U)] <- inf*sign(x_U[is.infinite(x_U)])
      
    # Linear Constraints:  
    A <- linCons[[1]]
    m <- nrow(A)
    b_L <- replicate(linCons[[2]], m)
    b_U <- replicate(linCons[[3]], m)
    b_L[is.infinite(b_L)] <- inf*sign(b_L[is.infinite(b_L)])
    b_U[is.infinite(b_U)] <- inf*sign(b_U[is.infinite(b_U)]) 
    
    # Optimize Portfolio:
    value <- amplLP(vec, x_L, x_U, A, b_L, b_U, control)
    
    # Return Value:
    value
}


# -------------------------------------------------------------------------


amplLP <- 
function(
    objective, x_L=NULL, x_U=NULL, A=NULL, b_L=NULL, b_U=NULL,
    control=list())
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Universal function wrapper for AMPL LP solvers
    
    # Arguments:
    #   objective - vec
    
    # FUNCTION:
    
    # Control List:
    ctrl <- amplLPControl()
    if (length(control) > 0)
        for (name in names(control)) ctrl[name] = control[name]
    control <- ctrl
    
    # Control Parameters:
    project <- control$project
    solver <- control$solver
    inf <- control$inf
    trace <- control$trace
    
    # Objective:
    c <- objective
    n <- length(vec)
    m <- nrow(A) 

    # Assign LP Model:
    .lpAssign(project, c, x_L, x_U, A, b_L, b_U , trace=FALSE)
    
    # Run AMPL:
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
    
    # Get Solver Message:
    Index <- grep("solve_message", out):length(out)
    message <- out[Index]
    
    # Version:
    version <- system(paste(solver, "-v"), intern=TRUE) 
    
    # Compose Results into a List:
    objective <- (c %*% solution)[[1, 1]]
    
    # Return Value:
    model <- capture.output(amplModelShow(project))
    run <- capture.output(amplModelShow(project))
    value = list(
        opt = list(solve=solve, model=model, run=run, out=out),
        solution = solution, 
        objective = objective,
        status = statusCode,
        message = statusMessage,
        solver = paste("AMPL", solver),
        version = version)
    class(value) <- c("solver", "list")
    value
}

    
# -------------------------------------------------------------------------


amplLPControl <- 
    function(solver="ipopt", project="ampl", inf=1e12, trace=FALSE) 
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns AMPL LP control parameter list
    
    # FUNCTION:
    
    # Return Value:
    list(solver=solver, project=project, inf=inf, trace=trace)
}


###############################################################################

