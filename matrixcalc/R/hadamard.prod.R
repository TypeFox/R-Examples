hadamard.prod <- function( x, y )
{
###
### this function calculates the Hadamard product of two matrices x and y.
### the matrices must the same row and column order
###
### Parameters
### x = a numeric matrix object
### y = a numeric matrix object
###
    if ( !is.numeric( x ) ) {
        stop( "argument x is not numeric" )
    }    
    if ( !is.numeric( y ) ) {
        stop( "argument y is not numeric" )
    }    
    if ( is.matrix( x ) ) {
        Xmat <- x
    }
    else {
        if ( is.vector( x ) ) {
            Xmat <- matrix( x, nrow=length(x), ncol=1 )
        }
        else {
            stop( "argument x is neither a matrix or a vector" )
        }
    }    
    if ( is.matrix( y ) ) {
        Ymat <- y
    }
    else {
        if ( is.vector( y ) ) {
            Ymat <- matrix( y, nrow=length(x), ncol=1 )
        }
        else {
            stop( "argument x is neither a matrix or a vector" )
        }
    }    
    if ( nrow( Xmat ) != nrow( Ymat ) )
        stop( "argumentx x and y do not have the same row order" )
    if ( ncol( Xmat ) != ncol( Ymat ) )
        stop( "arguments x and y do not have the same column order" )
    return( Xmat * Ymat )
}
