.onLoad <- function( libname, pkgname ) {
  library.dynam( "quadprog", package = "quadprog", lib.loc = .libPaths() )
    ## Solver plugin name (based on package name)
    if( ! pkgname %in% ROI_registered_solvers() ){
        ## Register solver methods here.
        ## One can assign several signatures a single solver method
        solver <- ROI:::get_solver_name( pkgname )
        ROI:::ROI_register_solver_method( signatures = ROI:::ROI_make_QP_signatures(),
                                          solver = solver,
                                          method =
            getFunction( "solve_QP", where = getNamespace(pkgname)) )
        ## Finally, for status code canonicalization add status codes to data base
        .add_status_codes()
    }
}

#.onUnload <- function( libpath ){
#    ROI:::ROI_deregister_solver_methods( solver = "glpk" )
#}
