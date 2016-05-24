rshortonly.test <- function( m, n=2, k=n, x.t=1, x.l=0, x.u=x.t, max.iter=1000 )
{
###
### This function generates m short only random portfolios with n investments where
### absolute weights are between given lower and upper bound and the sum of
### the absolute weights equals a given total sum.  Results are return as a matrix
###
### Arguments
### m = a positive integer value for the number of portfolios to be generated
### m = a positive integer value for the number of investments in a portfolio
### k = a positive integer value for the number of non-zero weights
### x.t = a numeric value for the sum of the investment weights
### x.l = a numeric value for the lower bound for each weight
### x.u = a numeric value for the upper bound for each weight
### max.iter = a positive integer value for the maximum iterations in the rejection method
###
### private function
###
    by.case <- function( case, number, size, total, lower, upper, iterations )
    {
        return( random.shortonly.test( n=number, k=size, x.t=total, x.l=lower, x.u=upper, max.iter=iterations ) )
    }
###
### validate arguments
###
    if ( x.l >= x.u ) {
        stop( "Argument x.l is greater than or equal to argument x.u" )
    }
    if ( x.l < 0 ) {
        stop( "Argument x.l is less than zero" )
    }
    if ( x.u < 0 ) {
        stop( "Argument x.u is less than zero" )
    }
    results <- lapply( 1:m, by.case, n, k, x.t, x.l, x.u, max.iter )
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
