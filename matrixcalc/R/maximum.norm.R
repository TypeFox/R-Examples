maximum.norm <- function( x )
{
###
### This function computes and returns the maximum norm of matrix x
###
### arguments
### x = a numeric matrix or vector
###
    if ( !is.numeric( x ) ) {
        stop( "argument x is not numeric" )
    }
    if ( is.matrix( x ) ) {
        Xmat <- x
    }
    else {
        if ( is.vector( x ) ) {
            X.mat <- x
        }
        else {
            stop( "argument is neither a matrix nor a vector" )
        }
    }
    norm <- max( abs( Xmat ) )
    return( norm )
}
