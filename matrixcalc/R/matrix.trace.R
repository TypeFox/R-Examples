matrix.trace <- function( x )
{
###
### this function returns the trace of the given square matrix
###
### parameters
### x = a numeric square matrix object
###
    if ( !is.square.matrix( x ) )
        stop( "argument x is not a square matrix" )
    return( sum( diag( x ) ) )
}
