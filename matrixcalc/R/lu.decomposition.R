lu.decomposition <- function( x )
{
###
### This function performs an LU decomposition of the given square matrix argument
### the results are returned in a list of named components.  The Doolittle decomposition
### method is used to obtain the lower and upper triangular matrices
###
### arguments
### x = a square numeric matrix
###
    if ( !is.square.matrix( x ) )
        stop( "argument x is not a square matrix" )
    if ( !is.numeric( x ) )
        stop( "argument x is not numeric" )
    n <- nrow( x )
    L <- matrix( 0, nrow=n, ncol=n )
    U <- matrix( 0, nrow=n, ncol=n )
    diag( L ) <- rep( 1, n )
    for ( i in 1:n ) {
        ip1 <- i + 1
        im1 <- i - 1
        for ( j in 1:n ) {
            U[i,j] <- x[i,j]
            if ( im1 > 0 ) {
                for ( k in 1:im1 ) {
                    U[i,j] <- U[i,j] - L[i,k] * U[k,j]
                }
            }
        }
        if ( ip1 <= n ) {
            for ( j in ip1:n ) {
                L[j,i] <- x[j,i]
                if ( im1 > 0 ) {
                    for ( k in 1:im1 ) {
                        L[j,i] <- L[j,i] - L[j,k] * U[k,i]
                    }
                }
                if ( U[i,i] == 0 )
                    stop( "argument x is a singular matrix" )
                L[j,i] <- L[j,i] / U[i,i]
            }    
        }
    }
    result <- list( L=L, U=U )
    return( result )
}
