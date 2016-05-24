rshortonly <- function( m, n=2, k=n, x.t=1, x.l=0, x.u=x.t, max.iter=1000 )
{
###
### This function generates m short only random portfolios with n investments where
### weights are between given lower and upper bound and the sum of
### the weights equals a given total sum.  Results are return as a matrix
###
### Arguments
### m = a positive integer value for the number of portfolios to be generated
### m = a positive integer value for the number of investments in a portfolio
### x.t = a numeric value for the sum of the absolute value of the investment weights
### x.l = a numeric value for the lower bound for the absolute value of each weight
### x.u = a numeric value for the upper bound for the absolute value of each weight
### max.iter = a positive integer value for the maximum iterations in the rejection method
###
### private function
###
    by.case <- function( case, number, size, total, lower, upper, iterations )
    {
        return( random.shortonly( n=number, k=size, x.t=total, x.l=lower, x.u=upper, max.iter=iterations ) )
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
    weights <- t( sapply( 1:m, by.case, n, k, x.t, x.l, x.u, max.iter ) )
    return( weights )
}
