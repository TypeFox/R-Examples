is.non.singular.matrix <- function( x, tol=1e-8 )
{
###
### This function determines if a square matrix is non singular by
### evaluating the determinant
###
### arguments
### x = a numeric square matrix
### tol = tolerance level for zero
###
    if (!is.square.matrix( x ) )
        stop( "argument x is not a square matrix" )
    if ( !is.numeric( x ) )
        stop( "argument x is not a numeric matrix" )
    det.x <- det( x )
    return( abs( det.x ) >= tol )
}
