fibonacci.matrix <- function( n )
{
###
### This function returns the order n + 1 Fibonacci matrix which
### is a square matrix derived from the Fibonacci sequence
###
### argument
### n = a positive integer
###
    if ( !is.numeric( n ) )
        stop( "argument n is not numeric" )
    if ( !is.vector( n ) )
        stop( "argument n is not the proper data type" )
    if ( length( n ) > 1 )
        stop( "argument n is not a scalr number" )
    if ( n != trunc( n ) )
        stop( "argument n is not an integer" )
    if ( n < 1 )
        stop( "argument n is less than 2" )
    np1 <- n + 1
    np2 <- n + 2
    f <- rep( 0, np2 )
    f[1] <- 1
    f[2] <- 1
    j <- 3
    while ( j <= np2 ) {
        f[j] <- f[j-1] + f[j-2]
        j <- j + 1
    }
    F <- matrix( 0, nrow=np1, ncol=np1 )
    for ( i in 0:n ) {
        ip1 <- i + 1
        for ( j in 0:n ) {
            jp1 <- j + 1
            if ( i - j + 1 >= 0 ) {
                F[ip1,jp1] <- f[i - j + 2]
            }
        }
    }
    return( F )
}
