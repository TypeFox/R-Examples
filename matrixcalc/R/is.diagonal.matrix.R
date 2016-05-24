is.diagonal.matrix <- function( x, tol=1e-8 )
{
###
### this function returns TRUE if the off diagonal elements in absolute
### value are less than the given tolerance and FALSE otherwise
###
### arguments
### x = a numeric square matrix
### tol = tolerance level for zero
###
    if (!is.square.matrix( x ) )
        stop( "argument x is not a square matrix" )
    if ( !is.numeric( x ) )
        stop( "argument x is not a numeric matrix" )
    y <- x
    diag( y ) <- rep( 0, nrow( y ) )
    return( all( abs( y ) < tol ) )
}
