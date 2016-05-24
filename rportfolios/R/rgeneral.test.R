rgeneral.test <- function(m,  n=2, k=n, p=0.5, x.u=1/k )
{
###
### This function generates m random general portfolios with n investments each
### Each portfolio has k positions that can be positive or negative.
### The probability that a given non-zero position is
### positive is p.  The maximum absolute exposure is x.u which has 1 / k as the default
###
### Arguments
### m = a positive integer value for the number of portfolios to be generated
### n = a positive integer value for the number of investments in a portfolio
### p = a positive numeric value for the probability that an investment has positive weight
### x.u = a positive numeric value for the maximum absolute exposure to an investment
###
    if ( n <= 0 ) {
        stop( "Argument n is not positive" )
    }    
    if ( ( p < 0 ) || ( p > 1 ) ) {
        stop( "Argument p is not a valid probability" )
    }
    if ( x.u <= 0 ) {
        stop( "Argument x.u is not positive" )
    }
    if ( m <= 0 ) {
        stop( "Argument m is not positive" )
    }
###
### private function
###
    by.case <- function( case, number, size, probability, limit )
    {
        return( random.general.test( n=number, k=size, p=probability, x.u=limit ) )
    }
    results <- lapply( 1:m, by.case, n, k, p, x.u )
###
### separate the investment weights and iterations into a matrix and vector
###
    xmatrix <- matrix( 0, nrow=m, ncol=n )
    iters <- rep( 0, m )
    for ( case in 1:m ) {
        thisResult <- results[[case]]
        iters[case] <- thisResult$iter
        xmatrix[case,] <- thisResult$x
    }
###
### create a new result list
###
    result <- list( xmatrix=xmatrix, iters=iters )
    return( result )
}
