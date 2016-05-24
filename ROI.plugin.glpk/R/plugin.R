## ROI plugin: GLPK
## based on Rglpk interface

## BASIC SOLVER METHOD
solve_OP <- function( x, control ){
    if(is.null(control))
        control <- list()
    ## ROI has its own translation table, thus we do not need to canonicalize in Rglpk
    control$canonicalize_status = FALSE
    if( all(ROI::types(x) == "C") )
        out <- .solve_LP( x, control )
    else
        out <- .solve_MILP( x, control )
    out
}

## SOLVER SUBMETHODS
.solve_LP <- function( x, control ) {
    solver <- ROI:::get_solver_name( getPackageName() )
    out <- Rglpk_solve_LP( terms(objective(x))[["L"]],
                           constraints(x)$L,
                           constraints(x)$dir,
                           constraints(x)$rhs,
                           bounds = bounds(x),
                           max = x$maximum,
                           control = control)
    ## FIXME: keep oiginal solution (return value)
    ROI:::canonicalize_solution( solution = out$solution,
                                 optimum = out$optimum,
                                 status = out$status,
                                 solver = solver )
}

.solve_MILP <- function( x, control ) {
    solver <- ROI:::get_solver_name( getPackageName() )
    out <- Rglpk_solve_LP( terms(objective(x))[["L"]],
                           constraints(x)$L,
                           constraints(x)$dir,
                           constraints(x)$rhs,
                           bounds = bounds(x),
                           types = types(x),
                           max = x$maximum,
                           control = control)
    ROI:::canonicalize_solution( solution = out$solution,
                                 optimum = out$optimum,
                                 status = out$status,
                                 solver = solver )
}

## STATUS CODES
.add_status_codes <- function(){
    ## GLPK
    ## from GLPK 4.34 reference manual and glpk.h (symbol, code, message)
    ## FIXME: change in solver interface, canonicalization now done in ROI
    solver <- ROI:::get_solver_name( getPackageName() )
    ROI:::add_status_code_to_db(solver,
                                1L,
                                "GLP_UNDEF",
                                "Solution is undefined."
                                )
    ROI:::add_status_code_to_db(solver,
                                2L,
                                "GLP_FEAS",
                                "Solution is feasible."
                                )
    ROI:::add_status_code_to_db(solver,
                                3L,
                                "GLP_INFEAS",
                                "Solution is infeasible."
                                )
    ROI:::add_status_code_to_db(solver,
                                4L,
                                "GLP_NOFEAS",
                                "No feasible solution exists."
                                )
    ROI:::add_status_code_to_db(solver,
                                5L,
                                "GLP_OPT",
                                "Solution is optimal.",
                                0L
                                )
    ROI:::add_status_code_to_db(solver,
                                6L,
                                "GLP_UNBND",
                                "Solution is unbounded."
                                )
    invisible(TRUE)
}
