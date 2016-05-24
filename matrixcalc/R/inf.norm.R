inf.norm <- function( x )
{
###
### This function returns the p = Inf norm wich is the maximum absolute row sum
### of the matrix x
###
### argument
### x = a numeric vector or matrix
###
    if ( !is.numeric( x ) ) {
        stop( "argument x is not numeric" )
    }
    if ( is.matrix( x ) ) {
        Xmat <- x
    }
    else {
        if ( is.vector( x ) ) {
            Xmat <- matrix( x, nrow=length(x), ncol=1 )
        }
        else {
            stop( "argument x is neither a matrix nor a vector" )
        }
    }
    norm <- max( apply( abs( Xmat ), 1, sum ) )
    return( norm )
}
