print.migpd <- function( x , ... ){
    cat( "\nA collection of ", dim( x$data )[ 2 ] , "generalized Pareto models.\n" )
    conv <- sapply( x$models , function( x ) x$convergence )
    if ( sum( conv ) == 0 ) cat( "\nAll models converged.\n" )
    else cat( "\nWARNING: Not all models converged.\n" )
    invisible()
}