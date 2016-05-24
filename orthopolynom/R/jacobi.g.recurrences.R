jacobi.g.recurrences <- function( n, p, q, normalized=FALSE )
{
###
### This function returns a data frame with n+1 rows and four columns
### containing the coefficients c, d, e and f of the recurrence relations
### for the order k Jacobi polynomial, Gk(p,q,x),
### and for orders k=0,1,...,n
###
### Parameters
### n = integer highest polynomial order
### p = first polynomial parameter
### q = second polynomial parameter
### normalized = boolean value.  If true, recurrences are for normalized polynomials
###
    if ( n < 0 )
        stop( "negative highest polynomial order" )
    if ( n != round( n ) )
        stop( "highest polynomial order is not integer" )
    if ( ( p - q ) <= -1 )
        stop( "p minus q less than or equal to -1" )
    if ( q <= 0 )
        stop( "q less than or equal to 0" )
    np1 <- n + 1
    r <- data.frame( matrix( nrow=np1, ncol=4 ) )
    names( r ) <- c( "c", "d", "e", "f" )
###
### special case where p = 0
###
    if ( abs( p ) < 1e-6 ) {
    
        if ( normalized ) {
            norms <- sqrt( jacobi.g.inner.products( n + 1, p, q ) )
###
###         j = 0
###
            c <- 1
            d <- - q
            e <- 1
            f <- 0
            
            rho.j <- sqrt( 2 / ( q * ( 1 - q ) ) )
            
            r[1,"c"] <- c
            r[1,"d"] <- d * rho.j
            r[1,"e"] <- e * rho.j
            r[1,"f"] <- 0
###
###         j = 1
###
            c <- 1
            d <- ( q - 2 ) / 3
            e <- 1
            f <- - q * ( q - 1 ) / 2
            
            rho.j <- 6 / sqrt(  ( 1 + q ) * ( 2 - q ) )
            rho.jm1 <- 6 * sqrt( 2 / ( q * ( 2 - q ) * ( 1 - q * q ) ) )
            
            r[2,"c"] <- c
            r[2,"d"] <- d * rho.j
            r[2,"e"] <- e * rho.j
            r[2,"f"] <- f * rho.jm1
###
###         j > 1
###
            j <- 2
            k <- 3

            while ( j <= n ) {
            
                c <- pochhammer( 2 * j - 2, 4 ) * ( 2 * j - 1 )
                d <- -( 2 * j * ( j ) + q * ( - 1 ) ) * pochhammer(2 * j - 2, 3 )
                e <- pochhammer( 2 * j - 2, 4 ) * ( 2 * j - 1 )
                f <- j * ( j + q - 1 ) * ( j - 1 ) * ( j - q ) * ( 2 * j + 1 )
                
                rho.j   <- norms[k]   / norms[k+1]
                rho.jm1 <- norms[k-1] / norms[k+1]
                
                r[k,"c"] <- c
                r[k,"d"] <- d * rho.j
                r[k,"e"] <- e * rho.j
                r[k,"f"] <- f * rho.jm1
                
                j <- j + 1
                k <- k + 1
                
            } # end while

            return( r )
            
        } # end if normalized
        
        else {
        
###
###         j = 0
###
            c <- 1
            d <- - q
            e <- 1
            f <- 0
            
            r[1,"c"] <- c
            r[1,"d"] <- d
            r[1,"e"] <- e
            r[1,"f"] <- 0
###
###         j = 1
###
            
            c <- 1
            d <- ( q - 2 ) / 3
            e <- 1
            f <- - q * ( q - 1 ) / 2

            r[2,"c"] <- c
            r[2,"d"] <- d
            r[2,"e"] <- e
            r[2,"f"] <- f
###
###         j > 1
###
            j <- 2
            k <- 3

            while ( j <= n ) {
            
                c <- pochhammer( 2 * j - 2, 4 ) * ( 2 * j - 1 )
                d <- -( 2 * j * ( j ) + q * ( - 1 ) ) * pochhammer(2 * j - 2, 3 )
                e <- pochhammer( 2 * j - 2, 4 ) * ( 2 * j - 1 )
                f <- j * ( j + q - 1 ) * ( j - 1 ) * ( j - q ) * ( 2 * j + 1 )
                
                r[k,"c"] <- c
                r[k,"d"] <- d
                r[k,"e"] <- e
                r[k,"f"] <- f
                
                j <- j + 1
                k <- k + 1
                
            } # end while
            
            return( r )
            
        } # end else
            
    } # end if abs( p )
    
###
### general case
###
    j <- 0
    k <- 1
    if ( normalized ) {
        norms <- sqrt( jacobi.g.inner.products( n + 1, p, q ) )
        while ( j <= n ) {
            if ( j == 0 ) {
                c <- 1
                d <- - q / ( p + 1 )
                e <- 1
                f <- 0
            }
            else {
                c <- pochhammer( 2 * j + p - 2, 4 ) * ( 2 * j + p - 1 )
                d <- -( 2 * j * ( j + p ) + q * ( p - 1 ) ) * pochhammer(2 * j + p - 2, 3 )
                e <- pochhammer( 2 * j + p - 2, 4 ) * ( 2 * j + p - 1 )
                f <- j * ( j + q - 1 ) * ( j + p - 1 ) * ( j + p - q ) * ( 2 * j + p + 1 )
            }
            
            rho.j <- norms[k] / norms[k+1]
            
            r[k,"c"] <- c
            r[k,"d"] <- d * rho.j
            r[k,"e"] <- e * rho.j
            
            if ( j == 0 )
                r[k,"f"] <- 0
            else {
                if ( k == 1 )
                    r[k,"f"] <- 0
                else {
                    rho.jm1 <- norms[k-1] / norms[k+1]
                    r[k,"f"] <- f * rho.jm1
                }    
            }
            j <- j + 1
            k <- k + 1
        }
        return( r )
    }
    else {
        while ( j <= n ) {
            if ( j == 0 ) {
                c <- 1
                d <- - q / ( p + 1 )
                e <- 1
                f <- 0
            }
            else {
                c <- pochhammer( 2 * j + p - 2, 4 ) * ( 2 * j + p - 1 )
                d <- -( 2 * j * ( j + p ) + q * ( p - 1 ) ) * pochhammer(2 * j + p - 2, 3 )
                e <- pochhammer( 2 * j + p - 2, 4 ) * ( 2 * j + p - 1 )
                f <- j * ( j + q - 1 ) * ( j + p - 1 ) * ( j + p - q ) * ( 2 * j + p + 1 )
            }
            r[k,"c"] <- c
            r[k,"d"] <- d
            r[k,"e"] <- e
            r[k,"f"] <- f
            j <- j + 1
            k <- k + 1
        }
        return( r )
    }
    return( NULL )
}
