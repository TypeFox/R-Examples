H.matrices <- function( r, c=r )
{
###
### This function returns a list of lists.  The length of the high
### level list is r.  Each component i of the high level list is
### a list of c components.  Each sub-component j of main component i
### is an r by c matrix.  Each matrix is the outer product of column
### i in the r by r identity and column j of the c by c identity matrix
###
### Arguments
### r = a positive integer for the number of rows in each matrix
### c = a positive integer for the number of columns in each matrix
###
    if ( missing( r ) )
        stop( "argument r is missing" )
    if ( !is.numeric( r ) )
        stop( "argument r is not numeric" )
    if ( r != trunc( r ) )
        stop( "argument r is not an integer" )
    if ( r < 2 )
        stop( "argument r is less than 2" )
    if ( !is.numeric( c ) )
        stop( "argument c is not numeric" )
    if ( c != trunc( c ) )
        stop( "argument c is not an integer" )
    if ( c < 2 )
        stop( "argument c is less than 2" )
    Ir <- diag( rep( 1, r ) )
    Ic <- diag( rep( 1, c ) )
    H <- list()
    for ( i in 1:r) {
        H[[i]] <- list()
        for ( j in 1:c ) {
            H[[i]][[j]] <- Ir[i,] %o% Ic[j,]
        }
    }
    return( H )
}
