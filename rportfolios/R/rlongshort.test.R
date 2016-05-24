rlongshort.test <- function( m, n=2, k=n, x.t.long=1, x.t.short=x.t.long , max.iter=2000, eps=1e-3 )
{
###
### This function generates m long short random portfolios with n investments where
### gross notional exposure is x.t.long + and x.t.short and the net notional
### exposure is x.t.long - x.t.short.  Results are returned as a matrix.
###
### Arguments
### m = a positive integer value for the number of portfolios to be generated
### m = a positive integer value for the number of investments in a portfolio
### k = a positive integer value for the number of non-zero positions
### x.t.long = a numeric value for the sum of the long exposures in the long portfolio
### x.t.short = a numeric value for the sum of the short exposures in the short portfolio
### max.iter = a positive integer value for the maximum iterations in the rejection method
### eps = a positive numeric value for the acceptance rejection method based on gross notional exposure
###
### private function
###
    by.case <- function( case, number, size, total.long, total.short, 
        iterations, epsilon )
    {
        return( random.longshort.test( n=number, k=size, x.t.long=total.long,
            x.t.short=total.short, 
            max.iter=iterations, eps=epsilon ) )
    }
###
### validate arguments
###
    if ( x.t.long < 0 ) {
        stop( "Argument x.t.long is less than zero" )
    }
    if ( x.t.short < 0 ) {
        stop( "Argument x.t.short is less than zero" )
    }
    results <- t( lapply( 1:m, by.case, n, k, x.t.long, x.t.short, max.iter, eps ) )
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
