fix.fp <- function( z, tol=0.0001 )
{
###
### this function rounds floating pt values that are within the specified tolerance
### of an integer
###
### Parameters
### z = a numeric floating point value
### tol = the numeric floating point tolerance relative to an integer
###
    x <- 1
    m <- round( -log10( tol ) ) - 1
    mm1 <- m - 1
    for ( n in 1:m ) {
        if ( n > 1 )
            x <- 10 * x
        f <- x * z
        i <- round( f )
        j <- i - mm1
        for ( k in 1:m ) {
            g <- 1.0 * ( j + k )
            if ( abs( f - g ) <= tol )
                return ( g / x )
        }
    }
    return( z )
}
