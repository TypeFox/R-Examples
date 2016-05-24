svd.inverse <- function( x )
{
###
### This function returns the Penrose-Moore pseudo inverse matrix
### for the given n by p matrix A using singular value decomposition
###
### arguments
### x = a numeric matrix or vector
### nu = positive integer for the number of left singular vectors to be computed
### nv = positive integer for the number of right singular vectors to be computed
### LINKPAC = logical value to determine if LINKPACK should be used.
###
    if ( !is.numeric( x ) ) {
        stop( "argument x is not numeric" )
    }
    if ( is.matrix( x ) ) {
        Xmat <- x
    }
    else {
        if ( is.vector( x ) ) {
            Xmat <- matrix( x, nrow=length( x ), ncol=1 )
        }
        else {
            stop( "argument x is neither a vector nor a matrix" )
        }
    }
    svdXmat <- svd( Xmat )
    U <- svdXmat$u
    d <- svdXmat$d
    if ( any( d == 0 ) ) {
        stop( "x has at least one zero singular value" )
    }
    if ( length( d ) == 1 ) {
        Dinv <- matrix( 1 / d, nrow=1 )
    }
    else {
        Dinv <- diag( 1 / d )
    }    
    V <- svdXmat$v
    Xinv <- V %*% Dinv %*% t( U )
    return( Xinv )
}
