requal <- function( m, n=2, k, x.t=1 )
{
###
### This function generates m long only random portfolios with k investments where
### there are n non-zero equal weights that sum to x.t
###
### Arguments
### m = a positive integer value for the number of portfolios to be generated
### n = a positive integer value for the number of investments in a portfolio
### k = a positive integer value for the number of non zero investments
###
### private function
###
    by.case <- function( case, number, size, total )
    {
        return( random.equal( n=number, k=size, x.t=total ) )
    }
###
### validate arguments
###
    weights <- t( sapply( 1:m, by.case, n, k, x.t ) )
    return( weights )
}
