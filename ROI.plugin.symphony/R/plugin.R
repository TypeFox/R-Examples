## ROI plugin: SYMPHONY
## based on Rsymphony interface

## BASIC SOLVER METHOD
solve_OP <- function( x, control ){
    if( all(ROI::types(x) == "C") )
        out <- .solve_LP( x, control )
    else
        out <- .solve_MILP( x, control )
    out
}

## SOLVER SUBMETHODS
.solve_LP <- function( x, control ) {
    solver <- ROI:::get_solver_name( getPackageName() )
    ## handle control args
    main_args <- list( obj = terms(objective(x))[["L"]],
                       mat = constraints(x)$L,
                       dir = constraints(x)$dir,
                       rhs = constraints(x)$rhs,
                       bounds = bounds(x),
                       max = x$maximum )
    ## handle STMs (shouldn't this be done on the Rsymphony side?)
    if( slam::is.simple_triplet_matrix(main_args$obj) )
        main_args$obj <- as.matrix( main_args$obj )
    foo_args <- names( as.list(args(Rsymphony_solve_LP)) )
    control_args <- as.list(control)[ names(control) %in% foo_args[!foo_args %in% c(names(main_args), "")] ]
    out <- do.call( Rsymphony_solve_LP, c(main_args, control_args) )
    ## FIXME: retranslate Rsymphony status code canonicalization
    status <- if( out$status == 0L )
        c(TM_OPTIMAL_SOLUTION_FOUND = 227L)
    else
        out$status
    ## FIXME: keep oiginal solution (return value)
    ROI:::canonicalize_solution( solution = out$solution,
                                 optimum = out$objval,
                                 status = status,
                                 solver = solver )
}

.solve_MILP <- function( x, control ) {
    solver <- ROI:::get_solver_name( getPackageName() )
    main_args <- list( obj = terms(objective(x))[["L"]],
                       mat = constraints(x)$L,
                       dir = constraints(x)$dir,
                       rhs = constraints(x)$rhs,
                       bounds = bounds(x),
                       types = types(x),
                       max = x$maximum )
    ## handle STMs (shouldn't this be done on the Rsymphony side?)
    if( slam::is.simple_triplet_matrix(main_args$obj) )
        main_args$obj <- as.matrix( main_args$obj )
    foo_args <- names( as.list(args(Rsymphony_solve_LP)) )
    control_args <- as.list(control)[ names(control) %in% foo_args[!foo_args %in% c(names(main_args), "")] ]
    out <- do.call( Rsymphony_solve_LP, c(main_args, control_args) )
    ## FIXME: retranslate Rsymphony status code canonicalization
    status <- if( out$status == 0L )
        c(TM_OPTIMAL_SOLUTION_FOUND = 227L)
    else
        out$status
    ROI:::canonicalize_solution( solution = out$solution,
                                 optimum = out$objval,
                                 status = status,
                                 solver = solver )
}

## STATUS CODES
.add_status_codes <- function( ) {

    ## SYMPHONY
    ## SYMPHONY 5.5.10 (no reference found yet)
    ## FIXME: better description of status in message

    solver <- ROI:::get_solver_name( getPackageName() )
    ROI:::add_status_code_to_db(solver,
                                0L,
                                "TM_OPTIMAL_SOLUTION_FOUND",
                                "(DEPRECATED) Solution is optimal. Compatibility status code will be removed in Rsymphony soon.",
                                0L
                                )
    ROI:::add_status_code_to_db(solver,
                                225L,
                                "TM_NO_PROBLEM",
                                "TM_NO_PROBLEM"
                                )
    ROI:::add_status_code_to_db(solver,
                                226L,
                                "TM_NO_SOLUTION",
                                "TM_NO_SOLUTION"
                                )
    ROI:::add_status_code_to_db(solver,
                                227L,
                                "TM_OPTIMAL_SOLUTION_FOUND",
                                "Solution is optimal.",
                                0L
                                )
    ROI:::add_status_code_to_db(solver,
                                228L,
                                "TM_TIME_LIMIT_EXCEEDED",
                                "TM_TIME_LIMIT_EXCEEDED"
                                )
    ROI:::add_status_code_to_db(solver,
                                229L,
                                "TM_NODE_LIMIT_EXCEEDED",
                                "TM_NODE_LIMIT_EXCEEDED"
                          )
    ROI:::add_status_code_to_db(solver,
                                230L,
                                "TM_TARGET_GAP_ACHIEVED",
                                "TM_TARGET_GAP_ACHIEVED"
                                )
    ROI:::add_status_code_to_db(solver,
                                231L,
                                "TM_FOUND_FIRST_FEASIBLE",
                                "TM_FOUND_FIRST_FEASIBLE"
                                )
    ROI:::add_status_code_to_db(solver,
                                232L,
                                "TM_FINISHED",
                                "TM_FINISHED"
                                )
    ROI:::add_status_code_to_db(solver,
                                233L,
                                "TM_UNFINISHED",
                                "TM_UNFINISHED"
                                )
    ROI:::add_status_code_to_db(solver,
                                240L,
                                "TM_FEASIBLE_SOLUTION_FOUND",
                                "TM_FEASIBLE_SOLUTION_FOUND"
                                )
    ROI:::add_status_code_to_db(solver,
                                235L,
                                "TM_SIGNAL_CAUGHT",
                                "TM_SIGNAL_CAUGHT"
                                )
    ROI:::add_status_code_to_db(solver,
                                -251L,
                                "TM_ERROR__NO_BRANCHING_CANDIDATE",
                                "TM_ERROR__NO_BRANCHING_CANDIDATE"
                                )
    ROI:::add_status_code_to_db(solver,
                                -252L,
                                "TM_ERROR__ILLEGAL_RETURN_CODE",
                                "TM_ERROR__ILLEGAL_RETURN_CODE"
                                )
    ROI:::add_status_code_to_db(solver,
                                -253L,
                                "TM_ERROR__NUMERICAL_INSTABILITY",
                                "TM_ERROR__NUMERICAL_INSTABILITY"
                                )
    ROI:::add_status_code_to_db(solver,
                                -254L,
                                "TM_ERROR__COMM_ERROR",
                                "TM_ERROR__COMM_ERROR"
                                )
    ROI:::add_status_code_to_db(solver,
                                -275L,
                                "TM_ERROR__USER",
                                "TM_ERROR__USER"
                                )
    invisible(TRUE)
}
