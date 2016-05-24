rgeneral <- function(m, n=2, k=n, p=0.5, x.u=1/k )
{
###
### This function generates m random general portfolios with n investments 
### Each portfolio has k positive or negative positions.
### The probability that a given non-zero position is positive is p.
### The maximum absolute exposure is x.u which has 1 / k as the default
###
### Arguments
### m = a positive integer value for the number of portfolios to be generated
### n = a positive integer value for the number of investments in each portfolio
### k = a positive integer value for the number of non-zero positions.
### p = a positive numeric value for the probability that a non-zero position
###     has positive weight
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
        return( random.general( n=number, k=size, p=probability, x.u=limit ) )
    }
    weights <- t( sapply( 1:m, by.case, n, k, p, x.u ) )
    return( weights )
}
