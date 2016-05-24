entrywise.norm <- function( x, p )
{
###
### This function computes the entrywise norm in which the absolute value of
### each element of matrix x is raised to the power p and the norm is the
### sum of these quantities
###
### arguments
### x = a numeric matrix or vector
### p = a real value for the power
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
            stop( "argument x is neither vector nor matrix" )
        }
    }
    if ( p == 0 ) {
        stop( "exponent p is zero" )
    }
    if ( is.infinite( p ) ) {
        return( maximum.norm( x ) )
    }    
    return( ( sum( abs( Xmat ) ^ p ) ) ^ ( 1 / p ) )
}
