slegendre.recurrences <- function( n, normalized=FALSE )
{
###
### This function returns a data frame with n+1 rows and four columns
### containing the coefficients c, d, e and f of the recurrence relations
### for the order k shifted Legendre polynomial, Pk(x), and orders k=0,1,...,n
###
### Parameter
### n = integer highest polynomial order
### normalized = a boolean value.  if true, the recurrences are for normalized polynomials
###
    if ( n < 0 )
        stop( "negative highest polynomial order" )
    if ( n != round( n ) )
        stop( "highest polynomial order is not an integer" )
    np1 <- n + 1
    r <- data.frame( matrix( nrow=np1, ncol=4 ) )
    names( r ) <- c( "c", "d", "e", "f" )
    j <- 0
    k <- 1
    if ( normalized ) {
        while ( j <= n ) {
            r[k,"c"] <-  j + 1
            r[k,"d"] <- - sqrt( ( 2 * j + 3 ) * ( 2 * j + 1 ) )
            r[k,"e"] <-  2 * sqrt( ( 2 * j + 3 ) * ( 2 * j + 1 ) )
            if ( j == 0 )
                r[k,"f"] <- 0
            else {
                r[k,"f"] <- j * sqrt( ( 2 * j + 3 ) / ( 2 * j - 1 ) )
            }
            j <- j + 1
            k <- k + 1
        }
        return( r )
    }
    else {
        while ( j <= n ) {
            r[k,"c"] <- j + 1
            r[k,"d"] <- -( 2 * j + 1 )
            r[k,"e"] <- 4 * j + 2
            r[k,"f"] <- j
            j <- j + 1
            k <- k + 1
        }
        return( r )
    }
    return( NULL )
}   
