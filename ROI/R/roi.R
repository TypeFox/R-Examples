################################################################################
## Package: ROI
## File:    roi.R
## Author:  Stefan Theussl
## Changed: 2011-10-03
################################################################################



################################################################################
## MAIN FUNCTION TO SOLVE OPTIMIZATION PROBLEMS USING ROI
################################################################################

##' This function makes a given solver (or searches for an appropriate
##' solver) solve the supplied optimization problem.
##'
##' @title Solve an Optimization Problem
##' @param x an optimization problem of class \code{"OP"}
##' @param solver a character vector specifying the solver to use.
##' @param control a list with additional control parameters for the
##' solver. This is solver specific so please consult the
##' corresponding documentation.
##' @param ... a list of control parameters (overruling those specified in \code{control}).
##' @return a list containing the solution and a message from the solver.
##' @author Stefan Theussl
##' @export
ROI_solve <- function( x, solver, control = NULL, ... ){

    dots <- list(...)
    control[names(dots)] <- dots

    x <- as.OP( x )

    ## handle the boundary case of no variables.
    ## FIXME: should also consider other
    #if( !length(terms(objective(x))$L) ) {
    #    y <- .solve_empty_OP(x)
    #    return(y)
    #}

    if(!length(objective(x)))
        stop("Cannot compute solution of empty objective function.")

    methods <- get_solver_methods( OP_signature(x) )
    SOLVE <- methods[[ solver ]]
    SOLVE( x, control )
}



################################################################################
## UTILITY FUNCTIONS TO QUERY SOLVERS
################################################################################

##' Retrieve all available solvers hosted in the given package repository.
##'
##' @title Solver Tools
##' @param ... arguments passed on to \code{\link{available.packages}}.
##' @return a named character vector.
##' @author Stefan Theussl
##' @export
ROI_available_solvers <- function( ... ){
    pkgs <- grep(.plugin_prefix(), rownames(utils::available.packages( ... )), value = TRUE)
    structure(pkgs, names = get_solver_name(pkgs))
}

##' Retrieve the names of all installed solvers.
##'
##' @title Solver Tools
##' @param ... arguments passed on to \code{\link{installed.packages}}.
##' @return a named character vector.
##' @author Stefan Theussl
##' @export
ROI_installed_solvers <- function( ... ){
    pkgs <- grep( .plugin_prefix(), rownames(utils::installed.packages( ... )),
                  value = TRUE )
    structure( pkgs, names = get_solver_name(pkgs) )
}

##' Retrieve the names of all registered solvers.
##'
##' @title Solver Tools
##' @return a named character vector.
##' @author Stefan Theussl
##' @export
ROI_registered_solvers <- function(){
    ## solvers registered
    get_solver_packages_from_db()
}



################################################################################
## HELPER FUNCTIONS (not exported)
################################################################################

## returns solver method
get_solver_methods <- function( signature ){
    entries <- do.call( solver_db$get_entries, as.list(signature) )
    solvers <- unlist(lapply( entries, function(x) x$solver ))
    structure( lapply(entries, function(x) x$FUN), names = solvers)
}

## returns available solvers from db
get_solvers_from_db <- function( ) {
    unique( solver_db$get_field_entries("solver", unlist = TRUE) )
}

## returns package names of available solvers from db
get_solver_packages_from_db <- function ( ){
    solvers <- get_solvers_from_db()
    structure( get_package_name(solvers), names = solvers )
}



################################################################################
## REGISTER NEW SOLVER METHODS
################################################################################

ROI_register_solver_method <- function( signatures, solver, method ){
    for( i in 1:nrow(signatures) )
        do.call(solver_db$set_entry, c(as.list(signatures[i,]),
                                             list(solver = solver),
                                             list(FUN = method)))

    invisible( TRUE )
}



################################################################################
## SIGNATURES
################################################################################

## make only one signature
ROI_make_signature <- function(...){
    dotargs <- list(...)
    required <- names(formals(OP))
    if( length(dotargs) < 2 )
        stop( sprintf("Signature element for '%s' and '%s' need to be given.",
                      required[1], required[2]) )
    length(dotargs) <- length(formals(OP))
    if( is.null(names(dotargs)) )
        names(dotargs) <- names(formals(OP))
    else {
        nam <- names(dotargs)
        nam[nam == ""] <-  names(formals(OP))[!(names(formals(OP)) %in% nam)]
        names(dotargs) <- nam
    }
    stopifnot( all(names(dotargs) %in% required) )

    ## FIXME: handle NULL case
    #.sort_types(unique(dotargs$types))
    #out <- dotargs[ names(dotargs)[-which(names(dotargs) == "types")] ]

    dotargs <- lapply(dotargs, function(x) if(is.null(x))
                                               FALSE
                                           else x)
    .make_signature(do.call(ROI_expand, dotargs))
}

ROI_expand <- function(...){
    base::expand.grid(..., stringsAsFactors = FALSE)
}

## make a set of signatures based on problem class
ROI_make_LP_signatures <- function()
    ROI_make_signature( objective = "L",
                        constraints = "L",
                        types = c("C"),
                        bounds = c(TRUE, FALSE),
                        maximum = c(TRUE, FALSE) )

ROI_make_QP_signatures <- function()
    ROI_make_signature( objective = "Q",
                        constraints = "L",
                        types = c("C"),
                        bounds = c(TRUE, FALSE),
                        maximum = c(TRUE, FALSE) )


ROI_make_MILP_signatures <- function()
    ROI_make_signature( objective = "L",
                        constraints = "L",
                        types = c("C", "I", "B", "CI", "CB", "IB", "CIB"),
                        bounds = c(TRUE, FALSE),
                        maximum = c(TRUE, FALSE) )


ROI_make_MIQP_signatures <- function()
    ROI_make_signature( objective = c("L", "Q"),
                        constraints = "L",
                        types = c("C", "I", "B", "CI", "CB", "IB", "CIB"),
                        bounds = c(TRUE, FALSE),
                        maximum = c(TRUE, FALSE) )

ROI_make_MIQCP_signatures <- function()
    ROI_make_signature( objective = c("L", "Q"),
                        constraints = c("L", "Q"),
                        types = c("C", "I", "B", "CI", "CB", "IB", "CIB"),
                        bounds = c(TRUE, FALSE),
                        maximum = c(TRUE, FALSE) )


.make_signature <- function( x ){
    stopifnot( ncol(x) == length(formals(OP)) )
    stopifnot( identical(colnames(x), names(formals(OP))) )
    types <- strsplit(as.character(x[["types"]]), "")
    types <- do.call(rbind, lapply( types, function(t) available_types() %in% t) )
    colnames(types) <- available_types()
    cbind(x[, colnames(x) != "types"], types)
}


.sort_types <- function(x){
    stopifnot( all(x %in% available_types()) )
    ord <- c(C = 1, I = 2, B = 3)
    ordered <- order(ord[x])
    x[ordered]
}
