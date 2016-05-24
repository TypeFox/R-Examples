is.idempotent.matrix <- function( x, tol=1e-8 )
{
###
### This function returns TRUE if the argument matrix is idempotent,
### that is, the product of the matrix with itself returns the same
### matrix within the prescribed tolerance and FALSE otherwise
###
### arguments
### x = a numeric square matrix
### tol = tolerance level for zero
###
    if (!is.square.matrix( x ) )
        stop( "argument x is not a square matrix" )
    if ( !is.numeric( x ) )
        stop( "argument x is not a numeric matrix" )
    return( all( abs( x %*% x - x ) < tol ) )
}
