random.longshort.test <- function( n=2, k=n, x.t.long=1, x.t.short=x.t.long, max.iter=2000, eps=1e-3 )
{
###
### This function constructs a random portfolio of weights in a long short portfolio.  It is
### algebraically a long only random portfolio plus a short only random portfolio. The gross
### notional exposure is x.t.long + x.t.short whereas the net notional exposure is
### x.t.long - x .t.short.
###
### Arguments
### n = a positive integer value for the number of investments
### k = a positive integer value for the number of non-zero positions
### x.t.long = a positive numeric value for the sum of the exposures in the long only leg
### x.t.short = a positive numeric value for the sum of the absolute exposures in the short only leg
### max.iter = a positive integer value for the number of iterations in the rejection method
### eps = a positive numeric value for the acceptance rejection method based on gross notional exposure
###
    iter <- 0
    more <- TRUE
    gross.t <- x.t.long + x.t.short
    while ( more ) {
        iter <- iter + 1
        if ( iter > 2 * max.iter ) {
            stop( "Maximum rejection iterations exceeded in random.longshort" )
        }    
        x.long <- random.longonly( n=k, k, x.t=x.t.long, x.l=0, x.u=x.t.long, max.iter )
        x.short <- random.shortonly( n=k, k, x.t=x.t.short, x.l=0, x.u=x.t.short, max.iter )
        x.longshort <- x.long + x.short
        gross.x <- sum( abs( x.longshort ) )
        if ( abs( gross.x - gross.t ) <= eps ) {
            investments <- sample( 1:n, k, replace=FALSE )
            x <- rep( 0, n )
            x[investments] <- x.longshort
            result <- list( x=x, iter=iter )
            return ( result )
        }
    }
    investments <- sample( 1:n, k, replace=FALSE )
    x <- rep( 0, n )
    x[investments] <- x.longshort
    result <- list( x=x, iter=iter )
    return( x.longshort )
}
