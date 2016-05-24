is.skew.symmetric.matrix <- function ( x, tol=1e-8 )
{
###
### this function returns TRUE if the matrix argument x is
### skew symmetric and FALSE otherwise
###
### argument
### x = a numeric square matrix
### tol = tolerance level for zero
###
    if ( !is.square.matrix( x ) )
        stop( "argument x is not a square matrix" )
    if ( !is.numeric( x ) )
        stop( "argument x is not a numeric matrix" )
    return( all( abs( x + t(x) ) < tol ) )
}
