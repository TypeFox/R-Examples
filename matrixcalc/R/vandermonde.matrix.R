vandermonde.matrix <- function( alpha, n )
{
###
### this function returns an m by n matrix of the powers of the alpha vector
###
### Parameters
### alpha = an m dimensional vector
### n = an integer
###
    if ( !is.vector( alpha ) )
        stop( "argument alpha is not a vector" )
    if ( !is.numeric( alpha ) )
        stop( "argument n is not a numeric vector" )
    m <- length( alpha )
    V <- matrix( 0, nrow=m, ncol=n )
    V[,1] <- rep( 1, m )
    j <- 2
    while ( j <= n ) {
       x <- alpha ^ ( j - 1 )
       V[,j] <- x
       j <- j + 1
   }
   return( V )
}
