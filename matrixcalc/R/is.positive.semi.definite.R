is.positive.semi.definite <- function( x, tol=1e-8 )
{
###
### this function determines if the given real symmetric matrix is positive semi definite
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
    if ( any( eigenvalues < 0 ) ) {
        return( FALSE )
    }
    return( TRUE )
}
