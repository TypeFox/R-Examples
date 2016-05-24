stirling.matrix <- function( n )
{
###
### this function constructs the Stirling matrix, a lower triangular matrix containing
### the Stirling numbers of the second kind
###
### argument
### n = a positive integer for the order of the matrix
###
    if ( !is.numeric( n ) )
        stop( "argument n is not numeric" )
    if ( !is.vector( n ) )
        stop( "argument n is not of the proper data type" )
    if ( length( n ) != 1 )
        stop( "argument n is not a scalar integer" )
    if ( n < 1 )
        stop( "arguemtn n is less than 1" )
    if ( n == 1 ) {
        return( S )
    }
    nm1 <- n - 1
    np1 <- n + 1
    S <- matrix( 0, nrow=np1, ncol=np1 )
    for ( i in 0:n ) {
        ii <- i + 1
        S[ii,ii] <- 1
    }
    for ( j in 0:n ) {
        jj <- j + 1
        S[1,jj] <- 0
        S[jj,1] <- 0
    }    
    for ( i in 1:nm1 ) {
        ii <- i + 1
        for ( j in 1:i ) {
            jj <- j + 1
            S[ii+1,jj] <- S[ii,jj-1] + j * S[ii,jj]
        }
    }
    S[1,1] <- 1
    return( S )
}
