ghermite.h.recurrences <- function( n, mu, normalized=FALSE )
{
###
### This function returns a data frame with n+1 rows and four columns
### containing the coefficients c, d, e and f of the recurrence relations
### for the order k generalized Hermite polynomial, Hk(mu,x), and orders k=0,1,...,n
###
### Parameter
### n = integer highest polynomial order
### mu = polynomial parameter
### normalized = a boolean value.  if true, the recurrences are for normalized polynomials
###
    if ( n < 0 )
        stop( "negative highest polynomial order" )
    if ( n != round( n ) )
        stop( "highest polynomial order is not an integer" )
    if ( mu <= -0.5 )
        stop( "parameter mu is less than or equal to -0.5" )
    if ( mu == 0 )
        return( hermite.h.recurrences( n, normalized ) )
    np1 <- n + 1
    r <- data.frame( matrix( nrow=np1, ncol=4 ) )
    names( r ) <- c( "c", "d", "e", "f" )
    j <- 0
    k <- 1
    two.mu <- 2 * mu
    if ( normalized ) {
        while ( j <= n ) {
###
###         if theta.j is zero, then j is even
###
            theta.j <- j - 2 * floor( j / 2)
###
###         j is even
###
            if ( theta.j == 0 ) {
                r[k,"c"] <- 1
                r[k,"d"] <- 0
                r[k,"e"] <- sqrt( 2 / ( j + two.mu + 1 ) )
                if ( j == 0 ) {
                    r[k,"f"] <- 0
                }
                else {
                    r[k,"f"] <- sqrt( j / ( j + two.mu + 1 ) )
                }
            }
###
###         j is odd
###
            else {
                r[k,"c"] <- 1
                r[k,"d"] <- 0
                r[k,"e"] <- sqrt( 2 / ( j + 1 ) )
                if ( j == 0 ) {
                    r[k,"f"] <- 0
                }
                else {
                    r[k,"f"] <- sqrt( ( j + two.mu ) / ( j + 1 ) )
                }
            }    
            j <- j + 1
            k <- k + 1
        }
        return( r )
    }
    else {
        while ( j <= n ) {
            theta.j <- j - 2 * floor( j / 2)
            r[k,"c"] <- 1
            r[k,"d"] <- 0
            r[k,"e"] <- 2
            r[k,"f"] <- 2 * ( j + two.mu * theta.j )
            j <- j + 1
            k <- k + 1
        }
        return( r )
    }
    return( NULL )
}   
