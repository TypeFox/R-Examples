frobenius.matrix <- function( n )
{
###
### this function returns an order n Frobenius matrix that is useful
### in numerical mathematics.
###
### argument
### n = a positive integer greater than 1
###
    if ( !is.numeric( n ) )
        stop( "argument n is not numeric" )
    if ( !is.vector( n ) )
        stop( "argument n is not the proper data type" )
    if ( length( n ) > 1 )
        stop( "argument n is not a scalr number" )
    if ( n != trunc( n ) )
        stop( "argument n is not an integer" )
    if ( n < 2 )
        stop( "argument n is less than 2" )
    F <- shift.down( diag( 1, n ) )
    nm1 <- n - 1
    j <- 0:nm1
    F[,n] <- (-1)^(nm1-j) * choose( n,j)
    return( F )
}
