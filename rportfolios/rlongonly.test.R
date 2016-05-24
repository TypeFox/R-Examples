rlongonly.test <- function( m, n=2, x.sum=1, x.lb=0, x.ub=x.sum, max.iter=1000 )
{
###
### This function generates a collection of random long only portfolios
### The weights are non-negative and sum to x.sum.  The maximum value
### That a weight can be is x.ub
###
### Arguments
### m = positive integer value for the number of portfolios in the collection
### n = positive integer value for the number of investments in a portfolio
### x.sum = a positive numerical value for the sum of the weights
### x.ub = a positive numerical value for the maximum investment weight
### max.iter = a positive integer value for the maximum number of iterations
###
    set.seed( 1234 )
###
### test the consistency of the number of investments, the maximum weight and
### the sum of the weights
###
    if ( ( n * x.ub ) <= x.sum ) {
        stop( "weight upper bound inconsistent with sum of weights and number of investments" )
    }
###
### run the simulation
###
    output <- lapply( 1:m, random.longonly.test, n, x.lb, x.sum, x.ub, max.iter )
    return( weights )
}
