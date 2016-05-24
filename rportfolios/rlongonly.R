rlongonly <- function( m, n=2, x.sum=1, x.lb=0, x.ub=x.sum, max.iter=1000 )
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
### private function to generate one portfolio
###
    this.longonly <- function( case, m, n, x.lb, x.ub, x.sum, max.iter )
    {
        print( paste( "case", case ) )
        print( paste( "x.lb", x.lb ) )
        print( paste( "x.ub", x.ub ) )
        print( paste( "x.sum", x.sum  ) )
        nm1 <- n - 1
        iter <- 0
        values <- rep( x.lb, n )
        more <- TRUE
        while( more ) {
            values[1] <- x.lb + ( x.ub - x.lb ) * runif( 1, min=0, max=1 )
            values[2:n] <- x.lb
            iter <- iter + 1
            if ( iter > max.iter ) {
                stop( "Maximum number of iterations exceeded" )
            }
            for ( i in 2:nm1 ) {
                im1 <- i - 1
                cumulative.x <- sum( values[1:im1] )
                upper.bound <- min( x.ub, x.sum - cumulative.x )
                if ( x.lb < upper.bound )
                    values[i] <- x.lb + ( upper.bound - x.lb ) * runif( 1, min=0, max=1 )
            }
            cumulative.x <- sum( values[1:nm1] )
            values[n] <- x.sum - cumulative.x
            if ( values[n] <= x.ub ) {
                indices <- sample( 1:n, n )
                weights <- values[indices]
                return( weights )
            }
        }
        cumulative.x <- cumsum( values[1:nm1] )
        values[n] <- x.sum - cumulative.x
        indices <- sample( 1:n, n )
        weights <- values[indices]
        return( weights )
    }
###
### test the consistency of the number of investments, the maximum weight and
### the sum of the weights
###
    if ( ( n * x.ub ) <= x.sum ) {
        stop( "weight upper bound inconsistent with sum of weights and number of investments" )
    }
###
### create the matrix for the random portfolio weights
###
    weights <- t( sapply( 1:m, this.longonly, m, n, x.lb, x.ub, x.sum, max.iter ) )
    return( weights )
}
