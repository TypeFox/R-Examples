is.square.matrix <- function( x )
{
###
### determines if the given matrix is a square matrix
###
### arguments
### x = a matrix object
###
    if ( !is.matrix( x ) )
        stop( "argument x is not a matrix" )
    return( nrow(x) == ncol(x) )
}
