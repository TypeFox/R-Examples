requal.test <- function( m, n=2, k, x.t=1 )
{
###
### This function generates m long only random portfolios with p investments where
### there are n non-zero equal weights that sum to one.
###
### Arguments
### m = a positive integer value for the number of portfolios to be generated
### n = a positive integer value for the number of investments in a portfolio
### k = a positive integer value for the number of non zero investments
### x.t = the sum of the investment weights
###
### private function
###
    by.case <- function( case, number, size, total )
    {
        return( random.equal.test( n=number, k=size, x.t=total ) )
    }
###
### validate arguments
###
    results <- lapply( 1:m, by.case, n, k, x.t )
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
