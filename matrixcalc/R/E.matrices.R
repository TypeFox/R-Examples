E.matrices <- function( n )
{
###
### This function creates a list of lists.  The number of components
### in the higher level list is n.  For each component i in the
### higher level list, the number of components in the sub-list
### is n.  Each component j of the sub-list is an n by n matrix
### e_i %*% t ( e_j ).  Each of the arrays is a vector in an
### n by n identity matrix
###
### argument
### n = a positive integer value greater than or equal to 2
###
    if( missing( n ) )
        stop( "argument n is missing" )
    if ( !is.numeric( n ) )
        stop( "argument n is not numeric" )
    if ( n != trunc( n ) )
        stop( "argument n is not an integer" )
    if ( n < 2 )
        stop( "argument n is less than 2" )
    I <- diag( rep( 1, n ) )
    E <- list()
    for ( i in 1:n ) {
        E[[i]] <- list()
        for ( j in 1:n ) {
            E[[i]][[j]] <- I[i,] %o% I[j,]
        }
    }
    return( E )
}
