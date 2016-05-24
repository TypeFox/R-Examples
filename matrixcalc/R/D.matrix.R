D.matrix <- function( n )
{
###
### This function constructs the linear transformation D that maps vech(A)
### to vec(A) when A is a symmetric matrix
###
### Arguments
### n = a positive integer for the order of a matrix
###
    if( missing( n ) )
        stop( "argument n is missing" )
    if ( !is.numeric( n ) )
        stop( "argument n is not numeric" )
    if ( n != trunc( n ) )
        stop( "argument n is not an integer" )
    if ( n < 2 )
        stop( "argument n is less than 2" )
    p <- n * ( n + 1 ) /2
    nsq <- n * n
    Dt <- matrix( 0, nrow=p, ncol=nsq )
    T <- T.matrices( n )
    u <- u.vectors( n )
    k <- u$k
    I <- u$I
    for ( j in 1:n ) {
        for ( i in j:n ) {
            Dt <- Dt + I[,k[i,j]] %*% t( vec( T[[i]][[j]] ) )
        }
    }
    return( t( Dt ) )
}
