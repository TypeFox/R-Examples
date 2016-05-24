jacobi.p.recurrences <- function( n, alpha, beta, normalized=FALSE )
{
###
### This function returns a data frame with n+1 rows and four columns
### containing the coeffieicnts c, d, e and f of the recurrence relations
### for the order k Jacobi polynomial, Pk(alpha,beta,x),
### and for orders k=0,1,...,n
###
### Parameters
### n = integer highest polynomial order
### alpha = first polynomial parameter
### beta = second polynomial parameter
### normalized = boolean value.  If true, recurrences are for normalized polynomials
###
    if ( ( abs( alpha ) < 1e-6 ) & ( abs( beta ) < 1e-6 ) )
        return( legendre.recurrences( n, normalized ) )
    if ( abs( alpha - beta ) < 1e-6 ) {
        alpha.prime <- alpha + 0.5
        return( gegenbauer.recurrences( n, alpha.prime, normalized ) )
    }    
    if ( n < 0 )
        stop( "negative highest polynomial order" )
    if ( n != round( n ) )
        stop( "highest polynomial order is not integer" )
    if ( alpha <= -1 )
        stop( "alpha less than or equal to -1" )
    if ( beta <= -1 )
        stop( "beta less than or equal to -1" )
    np1 <- n + 1
    r <- data.frame( matrix( nrow=np1, ncol=4 ) )
    names( r ) <- c( "c", "d", "e", "f" )
    j <- 0
    k <- 1
    ab <- alpha + beta
    ab.zero <- abs( alpha + beta ) < 1e-6
    aabb <- alpha * alpha - beta * beta
    if ( normalized ) {
        norms <- sqrt( jacobi.p.inner.products( n+1, alpha, beta ) )
        while ( j <= n ) {
            c <- 2 * ( j + 1 ) * ( j + ab + 1 ) * ( 2 * j + ab )
            d <- ( 2 * j + ab + 1 ) * aabb
            e <- pochhammer( 2 * j + ab, 3 )
            f <- 2 * ( j + alpha ) * ( j + beta ) * ( 2 * j + ab + 2 )
            if ( ab.zero && j == 0 ) {
                c <- 1
                d <- alpha
                e <- 1
                f <- 0
            }    
            if ( j == 0 ) {
                rho.j <- sqrt( ( ab + 3 ) / ( ( alpha + 1 ) * ( beta + 1 ) ) )
            }
            else {
                num <- ( j + 1 ) * ( 2 * j + ab + 3 ) * ( j + ab + 1 )
                den <- ( j + alpha + 1 ) * ( j + beta + 1 ) * ( 2 * j + ab + 1 )
                rho.j <- sqrt( num / den )
            }
#            rho.j <- norms[k] / norms[k+1]
            r[k,"c"] <- c
            r[k,"d"] <- d * rho.j
            r[k,"e"] <- e * rho.j
            if ( j == 0 ) {
                rho.jm1 <- 0
            }
            else {
                if ( j == 1 ) {
                    num <- 2 * ( ab + 5 ) * ( ab + 2 )
                    den <- (alpha + 1 ) * ( alpha + 2 ) * ( beta + 1 ) * ( beta + 2 )
                    rho.jm1 <- sqrt( num / den )
                }    
                else {
                    num <- j * ( j + 1 ) * ( 2 * j + ab + 3 ) * gamma( j + ab + 2 )
                    den <- ( 2 * j + ab - 1 ) * gamma( j + ab ) * 
                           ( j + alpha + 1 ) * ( j + beta + 1 ) * 
                           ( j + alpha ) * ( j + beta )
                    rho.jm1 <- sqrt( num / den )
                }
            }    
#            rho.jm1 <- norms[k-1] / norms[k+1]
            r[k,"f"] <- f * rho.jm1
            j <- j + 1
            k <- k + 1
        }
        return( r )
    }
    else {
        while ( j <= n ) {
            c <- 2 * ( j + 1 ) * ( j + ab + 1 ) * ( 2 * j + ab )
            d <- ( 2 * j + ab + 1 ) * aabb
            e <- pochhammer( 2 * j + ab, 3 )
            f <- 2 * ( j + alpha ) * ( j + beta ) * ( 2 * j + ab + 2 )
            if ( ab.zero  && j == 0) {
                c <- 1
                d <- alpha
                e <- 1
                f <- 0
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
