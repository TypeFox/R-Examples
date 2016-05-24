random.equal.test <- function( n=2, k, x.t=1 )
{
###
### This function returns an n by 1 vector of investment weights where
### there are only k non-zero weights each of which is x.t / k
###
### Arguments
### n = a positive integer value for the number of investments in the portfolio
### k = a positive integer value for the number of non-zero investments
### x.t = the sum of the weights
###
    if ( k >= n )
       stop( "Argument k is greater than or equal to n" )
    weight <- x.t / k
    investments <- sample( 1:n, k, replace=FALSE )
    x <- rep( 0, n )
    x[investments] <- weight
    result <- list( x=weight, iter=1 )
    return( result )
}
