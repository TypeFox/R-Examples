is.symmetric.matrix <- function( x )
{
###
### this function determines if the matrix is symmetric
###
### argument
### x = a numeric matrix object
###
    if ( !is.matrix( x ) ) {
        stop( "argument x is not a matrix" )
    }
    if ( !is.numeric( x ) ) {
        stop( "argument x is not a numeric matrix" )
    }    
    if ( !is.square.matrix( x ) )
        stop( "argument x is not a square numeric matrix" )
    return( sum( x == t(x) ) == ( nrow(x) ^ 2 ) )
}
