
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
#  rneosQP                      Rmetrics Interface for AMPL/NEOS QP solvers
#  neosQP                       Convenience wrapper for AMPL/NEOS QP solvers
#  neosQPControl                NEOS QP ontrol parameter list
###############################################################################
    

rneosQP <- 
    function(objective, lower=0, upper=1, linCons, control=list())
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Rmetrics Interface for AMPL/NEOS QP solvers
    
    # Arguments:
    #   objective - list(dvec=NULL, Dmat = NULL)
    
    # FUNCTION:
    
    # Control List:
    ctrl <- neosQPControl()
    if (length(control) > 0)
        for (name in names(control)) ctrl[name] = control[name]
    control <- ctrl
    
    # Control Parameters:
    solver <- control$solver
    category <- control$category
    project <- control$project
    inf <- control$inf
    trace <- control$trace    
    
    # General Settings:
    dvec <- objective$dvec
    Dmat <- objective$Dmat
    obj <- rbind(dvec, Dmat)
    
    # Box Constraints:
    replicate <- function(x, n) if(length(x) == 1) rep(x, n) else x
    n <- ncol(obj)
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
    value <- neosQP(objective, x_L, x_U, A, b_L, b_U, control)
    
    # Return Value:
    value
}


###############################################################################


neosQP <- 
    function(
    objective=list(dvec=NULL, Dmat=NULL), 
    x_L=NULL, x_U=NULL, A=NULL, b_L=NULL, b_U=NULL,
    control=list(), ...)
{   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Convenience wrapper for AMPL/NEOS QP solvers
    
    # Arguments:
    #   objective - list(dvec=NULL, Dmat = NULL)
    
    # FUNCTION:
    
    # Control List:
    ctrl <- neosQPControl()
    if (length(control) > 0)
        for (name in names(control)) ctrl[name] = control[name]
    control <- ctrl
    
    # Control Parameters:
    solver <- control$solver
    category <- control$category
    project <- control$project
    inf <- control$inf
    trace <- control$trace
    
    # Solver Settings:
    dvec <- objective$dvec
    Dmat <- objective$Dmat
    obj <- rbind(dvec, Dmat)
    n <- ncol(obj)
    m <- nrow(A)

    # Write AMPL Model File:
    amplModelOpen(project)
    model <- c(
        "param n ;",
        "param m ;",
        "param c{1..n} ;",
        "param F{1..n, 1..n} ;",        
        "param x_L{1..n} ;",
        "param x_U{1..n} ;", 
        "param A{1..m, 1..n} ;",
        "param b_L{1..m} ;",
        "param b_U{1..m} ;",   
        "var x{1..n};",
        "minimize Risk: sum {i in 1..n} x[i]*c[i] + 0.5*sum {i in 1..n} sum{j in 1..n} x[i]*F[i,j]*x[j] ;",
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
    amplDataAddVector(data="c", vector=dvec, project)
    amplDataAddMatrix(data="F", matrix=Dmat, project) 
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
    model <- paste(readLines(
      paste(project, "mod", sep=".")), sep = " ", collapse ="\n")
    data <- paste(readLines(
      paste(project, "dat", sep=".")), sep = " ", collapse ="\n")
    run <- paste(readLines(
      paste(project, "run", sep=".")), sep = " ", collapse ="\n")
    
    # Setup NEOS and AMPL Specifications:
    amplSpec <- list(model=model, data=data, commands=run, comments="NEOS")
    solverTemplate <- rneos::NgetSolverTemplate(
      category=category, solvername=solver, inputMethod="AMPL")
    xmls <- rneos::CreateXmlString(
      neosxml=solverTemplate, cdatalist=amplSpec)
    
    # Submit and Fetch NEOS Job:
    submittedJob <- rneos::NsubmitJob(
      xmlstring=xmls, user="rneos", interface="", id=0)
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
    objval <- (dvec %*% solution + 0.5 * solution %*% Dmat %*% solution)[[1, 1]]
    
    # Return Value:
    model <- capture.output(amplModelShow(project))
    run <- capture.output(amplRunShow(project))
    value = list(
        opt = list(solve=ans, model=model, run=run, out=out),
        solution = solution, 
        objective = objval,
        status = statusCode,
        message = statusMessage,
        solver = paste("AMPL", solver),
        version = version)
    class(value) <- c("solver", "list")
    value
}


# -----------------------------------------------------------------------------


neosQPControl <- 
    function(solver="ipopt", category="nco", project="neos", 
        inf=1e12, trace=FALSE) 
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Returns control parameter list
    
    # FUNCTION:
    
    # Control Parameter:
    control <- list(
        solver=solver, category=category, project=project, inf=inf, trace=trace)
        
    # Return Value:
    control
}



###############################################################################

