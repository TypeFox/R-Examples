
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
#  ramplNLP                 Rmetrics Interface for AMPL LP solvers 
#  amplNLP                  Convenience wrapper for AMPL LP solvers
#  amplControl              AMPL LP control parameter list         
################################################################################


ramplNLP <- 
    function(
        start, objective, lower=0, upper=1, amplCons, control=list(), ...)
{
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION:
    
    # Control List:
    ctrl = amplNLPControl()
    if (length(control) > 0)
        for (name in names(control)) ctrl[name] = control[name]
    control = ctrl
    
    # General Settings:
    if (length(start) == 1) {
        n = start
    } else {
        n = length(start)
    }
    
    # Box Constraints:
    if(length(lower) == 1) {
        par.lower = rep(lower, n)
    } else {
        par.lower = lower
    }
    if(length(upper) == 1) {
        par.upper = rep(upper, n)
    } else {
        par.upper = upper
    }
       
    # Controls:
    solver = control$solver
    project = control$project
    trace = control$trace
    
    # Add dots ...
    args <- list(...)
    m1 = names(match.call(expand.dots = TRUE))
    m2 = names(match.call(expand.dots = FALSE))
    Names = m1[!(m1 %in% m2)]
    amplDataOpen(project)
    amplDataAdd("n", data = n, type = "value", project)
    amplModelOpen(project)
    amplModelAdd("param n ;", project)
    count = 0
    for (a in args) {
        count = count + 1
        if(length(a) == 1) {
            type = "value"
        } else {
            if(is.matrix(a)) {
                type = "matrix"
            } else {
                type = "vector"
            }
        }
        if (type == "value") 
            amplModelAdd(paste("param ", Names[count], ";", sep = ""), project)
        if (type == "vector") 
            amplModelAdd(paste("param ", Names[count], "{1 ..", length(a), "};", sep = ""), project)
        if (type == "matrix") 
            amplModelAdd(paste("param ", Names[count], 
                "{1..", NCOL(a), ", 1..", NROW(a), "}", ";", sep=""), project)
        amplDataAdd(Names[count], data = a, type = type, project)
    }
    amplModelAdd(paste("param lower {1..", length(par.lower), "};", sep = ""), project) 
    amplModelAdd(paste("param upper {1..", length(par.upper), "};", sep = ""), project) 
    amplDataAdd("lower", data = par.lower, type = "vector", project)
    amplDataAdd("upper", data = par.upper, type = "vector", project)
    if(trace) amplDataShow(project)
    
    # Add Objective:
    amplModelAdd("var x{1..n} ;", project)
    amplModelAdd(paste( "minimize Function: ", objective, sep = ""), project)
    amplModelAdd(amplCons, project)
    if(trace) amplModelShow(project)  
    
    # Write Run File:
    solver = "ipopt"
    amplRunOpen(project)
    run <- c(
        "reset ;",
        paste("option solver ", solver, " ;", sep = ""),
        paste("model ", project, ".mod ;", sep = ""),
        paste("data ", project, ".dat ;", sep = ""),
        "solve ;",
        paste("display x > ", project, ".txt ;", sep = ""),
        "exit ;")
    amplRunAdd(run, project)
    if(trace) amplRunShow(project)
    
    # AMPL:
    command = paste("ampl", paste(project, "run", sep="."))
    solve = system(command, show.output.on.console = TRUE)

    # Read Result:
    file <- paste(project, "txt", sep = ".")
    ans = scan(file, what = character())
    ans = matrix(as.numeric(ans[-c(1:3, length(ans))]), byrow=TRUE, ncol = 2)
    index = sort(ans[,1], index.return = TRUE)$ix
    solution = ans[index, 2]
   
    # Return Value:
    value <- list(
        opt = args, 
        solution = solution, 
        objective = NA,
        status = NA,
        message = "none",
        solver = "amplNLP")
    class(value) = c("solver", "list")
    value      
}


# -----------------------------------------------------------------------------


amplNLP <-
    function()
{
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION:

    NA
}


###############################################################################


amplNLPControl <- 
    function(solver = "minos", project = "ampl", trace = FALSE)
{
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION:
    
    # Control Parameters:
    control <- list(
        solver = solver, 
        project = project, 
        trace = trace)
        
    # Return Value:
    control
}


###############################################################################
        
  