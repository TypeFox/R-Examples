
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
#  rneoslLP                     Rmetrics Interface for AMPL/NEOS LP solvers
#  neoslLP                      Convenience wrapper for AMPL/NEOS LP solvers
#  neoslLPControl               AMPL/NEOS LP control parameter list
###############################################################################


rneosLP <- 
    function(objective, lower=0, upper=1, linCons, control=list())
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Implements AMPL LP Interface
    
    # Arguments:
    #   objective - vec
    
    # FUNCTION:
    
    # Control List:
    ctrl <- neosLPControl()
    if (length(control) > 0) 
        for (name in names(control)) ctrl[name] = control[name]
    control <- ctrl
    
    # Controls:
    solver <- control$solver
    category <- control$category
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
    value <- neosLP(vec, x_L, x_U, A, b_L, b_U, control)
    
    # Return Value:
    value
}



# -------------------------------------------------------------------------


neosLP <- 
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
    ctrl <- neosLPControl()
    if (length(control) > 0)
        for (name in names(control)) ctrl[name] = control[name]
    control <- ctrl
    
    # Control Parameters:
    solver <- control$solver
    category <- control$category
    project <- control$project
    inf <- control$inf
    trace <- control$trace
    
    # Objective:
    c <- objective
    n <- length(vec)
    m <- nrow(A) 

    # Write AMPL Model File:
    amplModelOpen(project)
    model <- c(
        "param n ;",
        "param m ;",
        "param c{1..n} ;",      
        "param x_L{1..n} ;",
        "param x_U{1..n} ;", 
        "param A{1..m, 1..n} ;",
        "param b_L{1..m} ;",
        "param b_U{1..m} ;",   
        "var x{1..n};",
        "minimize Objective: sum {i in 1..n} x[i]*c[i] ;",
        "s.t. lower {i in 1..n}: x[i] >= x_L[i] ;",
        "s.t. upper {i in 1..n}: x[i] <= x_U[i] ;",   
        "s.t. linLower {j in 1..m}: sum{i in 1..n} A[j, i]*x[i] >= b_L[j] ;", 
        "s.t. linUpper {j in 1..m}: sum{i in 1..n} A[j, i]*x[i] <= b_U[j] ;", 
        NULL)
    amplModelAdd(model, project)
    if (trace) amplModelShow(project)
        
    # Write AMPL Data File:
    amplDataOpen(project)
    amplDataAddValue (data="n", value=n, project)
    amplDataAddValue (data="m", value=m, project)
    amplDataAddVector(data="c", vector=c, project)
    amplDataAddVector(data="x_L", vector=x_L, project)
    amplDataAddVector(data="x_U", vector=x_U, project)    
    amplDataAddMatrix(data="A", matrix=A, project)
    amplDataAddVector(data="b_L", vector=b_L, project)
    amplDataAddVector(data="b_U", vector=b_U, project)  
    if (trace) amplDataShow(project)
       
    # Write AMPL/NEOS RUN File:
    amplRunOpen(project)
    run <- c(
        "solve ;",
        "display x;",
        "display solve_result_num;", 
        "display solve_result;",
        "display solve_message;",
        "exit ;")
    amplRunAdd(run, project)
    if (trace) amplRunShow(project)
    
    # Get AMPL Files:
    model <- paste(readLines("ampl.mod"), sep = " ", collapse ="\n")
    data <- paste(readLines("ampl.dat"), sep = " ", collapse ="\n")
    run <- paste(readLines("ampl.run"), sep = " ", collapse ="\n")
    
    # Setup NEOS and AMPL Specifications:
    amplSpec <- list(model=model, data=data, commands=run, comments="NEOS")
    solverTemplate <- rneos::NgetSolverTemplate(category=category, solvername=solver, inputMethod="AMPL")
    xmls <- rneos::CreateXmlString(neosxml=solverTemplate, cdatalist=amplSpec)
    
    # Submit and Fetch NEOS Job:
    submittedJob <- rneos::NsubmitJob(xmlstring=xmls, user="rneos", interface="", id=0)
    ans <- rneos::NgetFinalResults(obj=submittedJob, convert=TRUE)
    out <- strsplit(ans@ans, split="\n")[[1]]
  
    # Get Weights:
    Index <- (grep("x .*. :=", out)+1):( grep("^;$", out)-1)
    Out <- out[Index]
    splits <- strsplit(paste(Out, collapse=" "), " ")[[1]]
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
    
    # Neos Job Version: 
    version <- out[1]
        
    # Compute Obective Function Value:
    objval <- (c %*% solution)[[1, 1]]
    
    # Return Value:
    model <- capture.output(amplModelShow(project))
    run <- capture.output(amplModelShow(project))
    value <- list(
        opt = list(solve=solve, model=model, run=run, out=out),
        solution = solution, 
        objective = objval,
        status = statusCode,
        message = statusMessage,
        solver = paste("AMPL", solver),
        version = version)
    class(value) <- c("solver", "list")
    value
}

    
# -------------------------------------------------------------------------


neosLPControl <- 
    function(solver="ipopt", category="lp", project="neos", 
    inf=1e12, trace=FALSE) 
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns control parameter list
    
    # FUNCTION:
    
    # Return Value:
    list(solver=solver, category=category, project=project, inf=inf, trace=trace)
}


###############################################################################

