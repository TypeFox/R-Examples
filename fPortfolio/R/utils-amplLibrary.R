
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
# FUNCTION:              DESCRIPTION:
#  .lpAssign              Assigns linear programming model
#  .qpAssign              Assigns quadratic programming model
################################################################################


.lpAssign <-
    function(project, c, x_L, x_U, A, b_L, b_U, solver="ipopt", trace=FALSE)
{
    # A function written by Diethelm Wuertz
    
    # Description:
    #   Assigns LP Model
    
    # Arguments:
    #   project - project name
    
    # FUNCTION:
    
    # Settings:
    n <- length(c)
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
       
    # Write AMPL RUN File:
    amplRunOpen(project)
    run <- c(
        paste("reset ;"),
        paste("option solver ", solver, " ;", sep = ""),
        paste("model ", project, ".mod ;", sep = ""),
        paste("data ", project, ".dat ;", sep = ""),
        paste("solve ;"),
        paste("display x > ", project, ".txt ;", sep = ""),
        paste("display solve_result_num > ", project, ".txt ;", sep = ""),
        paste("display solve_result > ", project, ".txt ;", sep = ""),
        paste("display solve_message > ", project, ".txt ;", sep = ""),
        paste("exit ;") )
    amplRunAdd(run, project)
    if (trace) amplRunShow(project)
    
    # Return Value:
    invisible()  
}


# -----------------------------------------------------------------------------


.qpAssign <-
    function(project, c, F, x_L, x_U, A, b_L, b_U, solver="ipopt", trace=FALSE)
{
    # A function written by Diethelm Wuertz
    
    # Description:
    #   Assigns LP Model
    
    # Arguments:
    #   project - project name
    
    # FUNCTION:
    
    # Settings:
    n <- length(c)
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
    amplDataAddVector(data="c", vector=c, project)
    amplDataAddMatrix(data="F", matrix=F, project) 
    amplDataAddVector(data="x_L", vector=x_L, project)
    amplDataAddVector(data="x_U", vector=x_U, project)    
    amplDataAddMatrix(data="A", matrix=A, project)
    amplDataAddVector(data="b_L", vector=b_L, project)
    amplDataAddVector(data="b_U", vector=b_U, project)  
    if (trace) amplDataShow(project)
       
    # Write AMPL RUN File:
    amplRunOpen(project)
    run <- c(
        paste("reset ;"),
        paste("option solver ", solver, " ;", sep = ""),
        paste("model ", project, ".mod ;", sep = ""),
        paste("data ", project, ".dat ;", sep = ""),
        paste("solve ;"),
        paste("display x > ", project, ".txt ;", sep = ""),
        paste("display solve_result_num > ", project, ".txt ;", sep = ""),
        paste("display solve_result > ", project, ".txt ;", sep = ""),
        paste("display solve_message > ", project, ".txt ;", sep = ""),
        paste("exit ;") )
    amplRunAdd(run, project)
    if (trace) amplRunShow(project)
    
    # Return Value:
    invisible()
}
    

################################################################################
    
       