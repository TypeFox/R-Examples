L.matrix <- function( n )
{
###
### This function constructs the elimination matrix as a mapping
### from vec(A) to vech(A)
###
### Arguments
### n = a positive integer value for the order of the matrix
###
    if ( missing( n ) )
        stop( "argument n is missing" )
    if ( !is.numeric( n ) )
        stop( "argument n is not numeric" )
    if ( n != trunc( n ) )
        stop( "argument n is not an integer" )
    if ( n < 2 )
        stop( "argument n is less than 2" )
    u <- u.vectors( n )
    E <- E.matrices( n )
    k <- u$k
    I <- u$I
    p <- n * ( n + 1 ) / 2
    nsq <- n * n
    L <- matrix( 0, nrow=p, ncol=nsq)
    for ( j in 1:n ) {
        for ( i in j:n ) {
            L <- L + I[,k[i,j]] %*% t( vec( E[[i]][[j]] ) )
        }
    }
    return( L )
}
