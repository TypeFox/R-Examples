spectral.norm <- function( x )
{
###
### This function computes the spectral norm for an m by n matrix x
### If x is a vector, then the L2 norm is returned
###
### arguments
### x = a numeric matrix or vector
###
    if ( !is.numeric( x ) ) {
        stop( "argument x is not numeric" )
    }
    if ( is.vector( x ) ) {
        return( sqrt( sum( x * x ) ) )
    }
    if ( !is.matrix( x ) ) {
        return( "argument x is not a matrix" )
    }
    A <- t(x) %*% x
    eigenA <- eigen( A )
    lambdaA <- eigenA$values
    maxLambdaA <- lambdaA[1]
    if ( maxLambdaA < 0 ) {
        stop( "t(x) %*% x is negative definite" )
    }
    return( sqrt( maxLambdaA ) )
}
