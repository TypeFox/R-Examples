direct.sum <- function( x, y )
{
###
### This function computes the direct sum of two vectors or
### matrices resulting in a block diagonal matrix
###
### Arguments
### x = a numeric matrix or vector
### y = a numeric matrix or vector
###
    matrices <- list()
    if ( is.vector( x ) ) {
        A <- matrix( x )
    }
    else {
        if ( is.matrix( x ) ) {
            A <- x
        }
        else {
            stop( "Argument x is not a matrix or vector" )
        }
    }
    if ( is.vector( y ) ) {
        B <- matrix( y )
    }
    else {
        if ( is.matrix( y ) ) {
            B <- y
        }
        else {
            stop( "Argument y is not a matrix or vector" )
        }
    }
    C <- rbind( cbind( A, matrix( 0, nrow=nrow(A), ncol=ncol(B) ) ),
                cbind( matrix( 0, nrow=nrow(B), ncol=ncol(A) ), B ) )
    return( C )
}
