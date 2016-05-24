random.longonly.test <- function( case, n=2, x.lb=0, x.sum=1, x.ub=x.sum, max.iter=1000 )
{
###
### This function generates one random long only portfolio for the given number of assets and other
### portfolio characteristics
###
### Arguments
### case = a positive integer value for the particular portfolio as a case
### n = a positive integer value which is the number of assets in the portfolio
### x.lb = a numeric value for the lower bound on the allocation to an asset
### x.sum = a numeric value for the sum of the allocations across all assets
### x.ub = a numeric value for the upper bound on the allocation to an asset
### max.iter = a positive integer for the maximum number of iterations in the rejection method loop
###
    nm1 <- n - 1
    iter <- 0
    values <- rep( x.lb, n )
    more <- TRUE
    while( more ) {
        U <- runif( n )
        values[1] <- x.lb + ( x.ub - x.lb ) * U[1]
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
                values[i] <- x.lb + ( upper.bound - x.lb ) * U[i]
        }
        cumulative.x <- sum( values[1:nm1] )
        values[n] <- x.sum - cumulative.x
        if ( values[n] <= x.ub ) {
            indices <- sample( 1:n, n )
            weights <- values[indices]
            result <- list( weights=weights, iter=iter )
            return( result )
        }
    }
    cumulative.x <- cumsum( values[1:nm1] )
    values[n] <- x.sum - cumulative.x
    indices <- sample( 1:n, n )
    weights <- values[indices]
    result <- list( weights=weights, iter=iter )
    return( result )
}
