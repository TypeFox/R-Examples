all.moments <- function( x, order.max=2, central=FALSE, absolute=FALSE, na.rm=FALSE )
{
    if( order.max < 2 )
        stop( "maximum order should be at least 2" )
    if ( is.matrix( x ) ) {
        n <- ncol( x )
        mu <- matrix( nrow=order.max+1, ncol=n )
        for ( order in 0:order.max )
            mu[order+1,] <- moment( x, order, central, absolute, na.rm )
        return( mu )    
    }
    else if ( is.vector( x ) ) {
        mu <- rep( 0, order.max+1 )
        for ( order in 0:order.max )
            mu[order+1] <- moment( x, order, central, absolute, na.rm )
        return( mu )
    }
    else if ( is.data.frame( x ) ) {
        n <- ncol( x )
        mu <- matrix( nrow=order.max+1, ncol=n )
        for ( order in 0:order.max )
            mu[order+1,] <- moment( x, order, central, absolute, na.rm )
        return( mu )
    }
    else
        return( all.moments( as.vector(x), order.max, central, absolute, na.rm ) )
}
