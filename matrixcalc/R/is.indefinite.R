is.indefinite <- function( x, tol=1e-8 )
{
###
### this function determines if the given real symmetric matrix is indefinite,
### that is, the matrix is neither positive definite nor negative definite.
###
### parameters
### x = a square numeric matrix object
### tol = tolerance level for zero
###
    if ( !is.square.matrix( x ) )
        stop( "argument x is not a square matrix" )
    if ( !is.symmetric.matrix( x ) )
        stop( "argument x is not a symmetric matrix" )
    if ( !is.numeric( x ) )
        stop( "argument x is not a numeric matrix" )
    eigenvalues <- eigen(x, only.values = TRUE)$values
    n <- nrow( x )
    for ( i in 1: n ) {
        if ( abs( eigenvalues[i] ) < tol ) {
            eigenvalues[i] <- 0
        }
    }
    if( any( eigenvalues > 0 ) && any( eigenvalues < 0 ) ) {
        return( TRUE )
    }
    return( FALSE )
}
