is.singular.matrix <- function( x, tol=1e-8 )
{
###
### this function returns TRUE if the matrix argument x
### is singular and FALSE otherwise
###
### argument
### x = a numeric square matrix
### tol = tolerance level for zero
###
    if (!is.square.matrix( x ) )
        stop( "argument x is not a square matrix" )
    if ( !is.numeric( x ) )
        stop( "argument x is not a numeric matrix" )
    det.x <- det( x )
    return( abs( det.x ) < tol )
}
