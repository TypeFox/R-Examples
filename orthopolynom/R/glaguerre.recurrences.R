glaguerre.recurrences <- function( n, alpha, normalized=FALSE )
{
###
### This function returns a data frame with n+1 and four columns
### containing the coefficients c, d, e and f of the recurrence relations
### for the order k generalized Laguerre polynomials, Lk(alpha,x),
### and for orders k=0,1,...n
###
### Parameters
### n = integer highest polynomial order
### alpha = polynomial parameter
### normalized = a boolean value.  If true, recurrences are for normalized polynomials
###
    if ( n < 0 )
        stop( "negative highest polynomial order" )
    if ( n != round( n ) )
        stop( "highest polynomial order is not integer" )
    if ( alpha <= -1 )
        stop( "alpha less than or equal to -1" )
    np1 <- n + 1
    r <- data.frame( matrix( nrow=np1, ncol=4 ) )
    names( r ) <- c( "c", "d", "e", "f" )
    j <- 0
    k <- 1
    if ( normalized ) {
        while ( j <= n ) {
            r[k,"c"] <- j + 1
            rho.j <- sqrt( ( j + 1 ) / ( alpha + j + 1 ) )
            r[k,"d"] <- ( 2 * j + alpha + 1 ) * rho.j
            r[k,"e"] <- - rho.j
            if ( j == 0 )
                r[k,"f"] <- 0
            else {
                r[k,"f"] <-  sqrt( j * ( j + 1 ) * ( j + alpha ) / (alpha + j + 1 ) )
            }
            j <- j + 1
            k <- k + 1
        }
        return( r )
    }
    else {
        while ( j <= n ) {
            r[k,"c"] <- j + 1
            r[k,"d"] <- 2 * j + alpha + 1
            r[k,"e"] <- -1
            r[k,"f"] <- j + alpha
            j <- j + 1
            k <- k + 1
        }
        return( r )
    }
    return( NULL )
}
