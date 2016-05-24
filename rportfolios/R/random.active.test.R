random.active.test <- function( x.b, x.g, k=length(x.b), max.iter=2000, eps=1e-3 )
{
###
### This function generates an active portfolio relative to the
### given benchmark portfolio.  It is the benchmark portfolio plus
### a notional neutral long short portfolio with the given gross
### notional amount.
###
### Arguments
### x.b = a numeric vector with the benchmark weights
### x.g = a numeric value for the gross notional amount of the long short portfolio
### max.iter = a positive integer value for the number of iterations in the rejection method
### eps = a positive numeric value for the acceptance rejection method based on gross notional exposure
###
    if ( !is.vector( x.b ) )
        stop( "Argument x.b is not a vector" )
    if ( !is.numeric( x.b ) )
        stop( "Argument x.b is not a numeric vector" )
    n <- length( x.b )
    if ( n == 1 )
        stop( "Argument x.b must be of length greater than 1" )
    if ( x.g <= 0 )
        stop( "Argument x.g is not positive" )
    result.ls <- random.longshort.test( n, k, x.t.long=x.g/2, x.t.short=x.g/2, max.iter, eps )    
    x.ls <- result.ls$x
    iter.ls <- result.ls$iter
    x.a <- x.b + x.ls
    result <- list( x=x.a, iter=iter.ls )
    return( result )
}
