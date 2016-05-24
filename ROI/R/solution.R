################################################################################
## Package: ROI
## File:    solution.R
## Author:  Stefan Theussl
## Changed: 2011-10-05
################################################################################



################################################################################
## Solution object
################################################################################

.solve_empty_OP <- function( x ) {
    ## Check whether constraints are satisfied (interpreting each lhs as
    ## empty sum with value 0):
    constraints <- split( constraints(x)$rhs, constraints(x)$dir )
    if( all(unlist(Map(function(dir, rhs) get(dir)(0, rhs),
                       names(constraints), constraints))) )
        make_OP_solution( double(), 0, 0L, solver = "ROI_NULL" )
    else
        make_OP_solution( double(), NA_real_, 2L, solver = "ROI_NULL" )
}

make_OP_solution <- function(solution, objval, status, solver, ...)
    structure( list(solution = solution,
                    objval = objval,
                    status = status),
              meta = list(solver = solver, ...),
              class = "OP_solution" )



################################################################################
## Methods on solution object
################################################################################

##' @noRd
##' @S3method print OP_solution
print.OP_solution <- function(x, ...){
    success <- x$status$code == 0
    if( !success ){
        writeLines( "No solution found." )
        writeLines( sprintf("The solver message was: %s", x$status$msg$message) )
    } else{
        writeLines( "Optimal solution found." )
    }
    writeLines( sprintf("The objective value is: %e", x$objval) )
}



################################################################################
## CANONICALIZER
################################################################################

canonicalize_status <- function( status, solver ){
    msg <- get_status_message_from_db( solver, status )
    list( code = msg$roi_code, msg = msg )
}

canonicalize_solution <- function( solution, optimum, status, solver, ... )
{
    status <- canonicalize_status( status, solver )
    make_OP_solution( solution = solution,
                      objval  = optimum,
                      status   = status,
                      solver   = solver, ... )
}
