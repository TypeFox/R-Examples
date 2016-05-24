schebyshev.u.recurrences <- function( n, normalized=FALSE )
{
###
### This function returns a data frame with n+1 rows and four columns
### containing the coefficients c, d, e and f of the recurrence relations
### for the order k shifted Chebyshev polynomial of the second kind, Uk(x),
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
    j <- 0
    k <- 1
    if ( normalized ) {
        while ( j <= n ) {
            r[k,"c"] <-   1
            r[k,"d"] <-  -2
            r[k,"e"] <-   4
            if ( j == 0 )
                r[k,"f"] <- 0
            else {
                r[k,"f"] <- 1
            }
            j <- j + 1
            k <- k + 1
        }
        return( r )
    }
    else {
        r$c <- rep(  1, np1 )
        r$d <- rep( -2, np1 )
        r$e <- rep(  4, np1 )
        r$f <- rep(  1, np1 )
        return( r )
    }
    return( NULL )
}
