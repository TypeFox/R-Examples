##' @noRd
##' @import registry methods
ROI_options <-
local({
    options <- list()
    function(option, value) {
        if (missing(option)) return(options)
        if (missing(value))
            options[[option]]
        else
            options[[option]] <<- value
    }
})

## STATUS_DB
## create registry object containing status codes
add_status_db_schema <- function( status_db ){
    status_db$set_field("solver",   type = "character", is_key = TRUE)
    status_db$set_field("code",     type = "integer",   is_key = TRUE)
    status_db$set_field("symbol",   type = "character")
    status_db$set_field("message",  type = "character")
    status_db$set_field("roi_code", type = "integer",   alternatives = 0:1)
    status_db
}
status_db <- registry( )
status_db <- add_status_db_schema( status_db )

## SOLVER_DB
add_solver_db_schema <- function( solver_db ){
    solver_db$set_field( "solver",      type = "character", is_key = TRUE )
    solver_db$set_field( "objective",   type = "character", validity_FUN = function(x) x %in% names(available_objective_classes()), is_key = TRUE)
    solver_db$set_field( "constraints", type = "character", validity_FUN = function(x) x %in% names(available_constraint_classes()), is_key = TRUE)
    for( type in available_types() )
        solver_db$set_field( type,      type = "logical", is_key = TRUE)
    solver_db$set_field( "bounds",      type = "logical", is_key = TRUE)
    solver_db$set_field( "maximum",     type = "logical", is_key = TRUE)
    solver_db$set_field( "FUN",         type = "function" )
    solver_db
}
solver_db <- registry( )
solver_db <- add_solver_db_schema( solver_db )

cross_validate_schema <- function( args, solver_db){
    stopifnot( all(args %in% c(solver_db$get_field_names(), "types")) )
    invisible( TRUE )
}

schema_valid <- cross_validate_schema( names(formals(OP)), solver_db )

.onLoad <- function( libname, pkgname ) {
    if( ! "ROI.plugin.nlminb" %in% ROI_registered_solvers() ){
        ## Register solver methods here.
        ## One can assign several signatures a single solver method
        ## DISABLED! see R code (solution of QP from examples.R in work not same as quadprog)
        #solver <- "nlminb"
        #ROI_register_solver_method( signatures = ROI_make_QP_signatures(),
        #                            solver = solver,
        #                            method =
        #                            getFunction( ".solve_QP_nlminb", where = getNamespace(pkgname)) )
        ## Finally, for status code canonicalization add status codes to data base
        #.add_nlminb_status_codes()
    }

    ## SET DEFAULTS: for the time being 'ROI_NULL' for solving empty
    ## OPs is the default solver
    ROI_options("default_solver", "ROI_NULL")
}

.onAttach <- function( libname, pkgname ){
    ## Search for all solvers in same library as ROI and register found solvers
    ## implicitely be running the corresponding .onLoad() function.
    solvers <- ROI_installed_solvers( lib.loc = libname )
    lapply( solvers, function( pkgname ){ nmspc <- tryCatch(getNamespace(pkgname), error = identity)
                                          if( !inherits(nmspc, "error") ){
                                              load <- methods::getFunction( ".onLoad", where = nmspc )
                                              load( libname = libname, pkgname = pkgname )
                                          }} )
    ## Startup messages
    packageStartupMessage( sprintf("%s: R Optimization Infrastructure", pkgname) )
    packageStartupMessage( sprintf("Registered solver plugins: %s.",
                                   paste(names(ROI_registered_solvers()), collapse = ", ")) )
    packageStartupMessage( sprintf("Default solver: %s.", ROI_options("default_solver")) )
}
