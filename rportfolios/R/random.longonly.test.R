random.longonly.test <- function(  n=2, k=n, x.t=1, x.l=0, x.u=x.t, max.iter=1000 )
{
###
### This function generates one random long only portfolio for the given number of assets and other
### portfolio characteristics
###
### Arguments
### n = a positive integer value which is the number of assets in the portfolio
### k = a positive integer value which is the number of non-zero weights
### x.t = a numeric value for the sum of the allocations across all assets
### x.l = a numeric value for the lower bound on the allocation to an asset
### x.u = a numeric value for the upper bound on the allocation to an asset
### max.iter = a positive integer for the maximum number of iterations in the rejection method loop
###
    if ( n < 2 ) {
        stop( "Argument n is less than 2 in random.longonly" )
    }
    if ( k > n ) {
        stop( "Argument k is greater than n" )
    }
###
### cardinality is less than n
###
    if ( k < n ) {
        thisResult <- random.longonly.test( n=k, k, x.t, x.l, x.u, max.iter )
        weights <- thisResult$x
        iter <- thisResult$iter
        investments <- sample( 1:n, k, replace=FALSE )
	x <- rep( 0, n )
	x[investments] <- weights
        result <- list( x=x, iter=iter )
        return( result )
    }    
    if ( x.t < 0 ) {
        stop( "Argument x.t is negative in random.longonly" )
    }
    if ( x.u < 0 ) {
        stop( "Argument x.u is negative in random.longonly" )
    }
    if ( x.l < 0 ) {
        stop( "Argument x.l is negative in random.longonly" )
    }    
    if ( n * x.u <= x.t ) {
        stop( "The product of argument n and x.u is less than or equal to argument x.t in random.longonly" )
    }
    if ( n * x.l >= x.t ) {
        stop( "The product of argument n and x.l is greater than or equal to argument x.t in random.longlonly" )
    }
    if ( x.l > 0 ) {
        total <- x.t - n * x.l
        upper <- x.u - x.l
        iterations <- max.iter
        thisResult <- random.longonly.test( n, x.t=total, x.l=0, x.u=upper, max.iter=iterations )
        x <- thisResult$x
        iter <- thisResult$iter
        x.p <- x + x.l
        result <- list( x=x.p, iter=iter )
        return( result)
    }    
    nm1 <- n - 1
    iter <- 0
    more <- TRUE
    while( more ) {
        values <- rep( 0, n )
        U <- runif( n )
        values[1] <- x.u * U[1]
        iter <- iter + 1
        for ( i in 2:nm1 ) {
            im1 <- i - 1
            cumulative.x <- sum( values[1:im1] )
            if ( cumulative.x < x.t ) {
                upper.bound <- min( x.u, x.t - cumulative.x )
                if ( upper.bound > 0 ) {
                    values[i] <- upper.bound * U[i]
                }    
            }    
        }
        cumulative.x <- sum( values[1:nm1] )
        values[n] <- x.t - cumulative.x
        if ( ( values[n] >= 0 ) && (values[n] <= x.u ) ) {
            indices <- sample( 1:n, n )
            weights <- values[indices]
            result <- list( x=weights, iter=iter )
            return( result )
        }
        if ( iter > max.iter ) {
            stop( "Maximum number of iterations exceeded in random.longonly" )
        }
    }
    cumulative.x <- sum( values[1:nm1] )
    values[n] <- x.t - cumulative.x
    indices <- sample( 1:n, n )
    weights <- values[indices]
    result <- list( x=weights, iter=iter )
    return( result )
}
