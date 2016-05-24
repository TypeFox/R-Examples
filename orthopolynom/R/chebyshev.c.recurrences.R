chebyshev.c.recurrences <- function( n, normalized=FALSE )
{
###
### This function returns a data frame with n+1 rows and four columns
### containing the coefficients c, d, e and f of the recurrence relations
### for the order k Chebyshev polynomial of the first kind, Ck(x),
### and for orders k=0,1,...,n
###
### Parameter
### n = integer highest order
### normalized = boolean value.  If true, recurrences are for normalized polynomials
###
    if ( n < 0 )
        stop( "negative highest polynomial order" )
    if ( n != round( n ) )
        stop( "highest polynomial order is not integer" )
    np1 <- n + 1
    r <- data.frame( matrix( nrow=np1, ncol=4 ) )
    names( r ) <- c( "c", "d", "e", "f" )
    if ( normalized ) {
        j <- 0
        k <- 1
        while ( j <= n ) {
            r[k,"c"] <- 1
            r[k,"d"] <- 0
            if ( j == 0 ) {
                r[k,"e"] <- 0.5 * sqrt( 2 )
                r[k,"f"] <- 0
            }
            else {
                r[k,"e"] <- 1
                if ( j == 1 ) {
                    r[k,"f"] <- sqrt( 2 )
                }
                else {
                    r[k,"f"] <- 1
                }    
            }
            j <- j + 1
            k <- k + 1
        }
        return( r )
    }
    else {
        r$c <- rep( 1, np1 )
        r$d <- rep( 0, np1 )
        r$e <- rep( 1, np1 )
        r$f <- rep( 1, np1 )
        r[1,"e"] <- 0.5
        return( r )
    }
    return( NULL )
}
