gegenbauer.recurrences <- function( n, alpha, normalized=FALSE )
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
    if ( alpha <= -0.5 )
        stop( "alpha less than or equal to -0.5" )
###
### alpha = 1.0
### special case is the Chebyshev polynomial of 
###the second kind U_k(x)
###
    if ( abs( alpha - 1 ) < 1e-6 )
        return( chebyshev.u.recurrences( n, normalized ) )
###
### alpha = 0.5
### special case is the Legendre polynomial P_k(x)
###
    if ( abs( alpha - 0.5 ) < 1e-6 )
        return( legendre.recurrences( n, normalized ) )
    np1 <- n + 1
    r <- data.frame( matrix( nrow=np1, ncol=4 ) )
    names( r ) <- c( "c", "d", "e", "f" )
    j <- 0
    k <- 1
###
### alpha = 0,0
### special case is related to the Chebyshev polynomial of 
### the first kind T_k(x)
###
    if ( abs( alpha ) < 1e-6 ) {
    
        if ( normalized ) {
        
            while( j <= n ) {
                r[k,"c"] <- j + 1
                r[k,"d"] <- 0
                if ( j == 0 ) {
                    r[k,"e"] <- sqrt( 2 )
                }
                else {
                    r[k,"e"] <- 2 * ( j + 1 )
                }
                if ( j == 0 ) {
                    r[k,"f"] <- 0
                }
                else {
                    if ( j == 1 ) {
                        r[k,"f"] <- 2 * sqrt( 2 )
                    }
                    else {
                        r[k,"f"] <- j + 1
                    }
                }    
                j <- j + 1
                k <- k + 1
                
            } # end while block    
        
        } # end if normalized block
        
        else {
        
            while ( j <= n ) {
                r[k,"c"] <- j + 1
                r[k,"d"] <- 0
                if ( j == 0 ) {
                    r[k,"e"] <- 2
                }
                else {
                    r[k,"e"] <- 2 * j
                }
                if ( j == 0 ) {
                    r[k,"f"] <- 0
                }
                else {
                    if ( j == 1 ) {
                        r[k,"f"] <- 2
                    }
                    else {
                        r[k,"f"] <- j - 1
                    }    
                }
                j <- j + 1
                k <- k + 1
                
            } # end while
            
        } # end else not normalized block
        
        return( r )
        
    } # end if block 
###
### general case
###
    two.alpha <- 2 * alpha
    
    if ( normalized ) {
    
        while ( j <= n ) {
            r[k,"c"] <- ( j + 1 )
            r[k,"d"] <- 0
            rho.j <- sqrt( ( alpha + j + 1 ) * ( j + 1 ) * gamma( two.alpha + j ) /
                             ( ( alpha + j ) * gamma( two.alpha + j + 1 ) ) )
            r[k,"e"] <- 2 * ( alpha + j ) * rho.j
            if ( j == 0 ) {
                rho.jm1 <- 0
            }    
            else {
                rho.jm1 <- sqrt( j * ( j + 1 ) * ( alpha + j + 1 ) * gamma( two.alpha + j - 1 ) /
                                   ( ( alpha + j - 1 ) * gamma( two.alpha + j + 1 ) ) )
            }
            r[k,"f"] <- ( j + 2 * alpha - 1 ) * rho.jm1
            j <- j + 1
            k <- k + 1
            
        } # end while block
        
    } # end if block
    
    else {
    
        while ( j <= n ) {
            r[k,"c"] <- j + 1
            r[k,"d"] <- 0
            r[k,"e"] <- 2 * ( j + alpha )
            r[k,"f"] <- j + 2 * alpha - 1
            j <- j + 1
            k <- k + 1
            
        } # end while block
        
    } # end of else block
    
    return( r )
}
