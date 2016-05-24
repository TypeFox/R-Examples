all.cumulants <- function(mu.raw)
{
    mu.central <- raw2central( mu.raw )
    if ( is.vector( mu.raw ) ) {
        order.max <- length( mu.raw ) - 1
        kappa <- rep( 0, order.max + 1 )
        kappa[2] <- mu.central[2]
        kappa[3] <- mu.central[3]
        n <- 3
        while( n <= order.max ) {
            np1 <- n + 1
            nm1 <- n - 1
            total <- mu.raw[np1]
            for ( k in 1:nm1 ) {
                km1 <- k - 1
                total <- total - choose(nm1,km1 ) * kappa[k+1] * mu.raw[n-k+1]
            }
            kappa[np1] <- total
            n <- n + 1
        }
        return( kappa )
    }
    else if ( is.matrix( mu.raw ) ) {
        order.max <- nrow( mu.raw ) - 1
        nvar <- ncol( mu.raw )
        kappa <- matrix( nrow=order.max+1,ncol=nvar )
        kappa[1,] <- rep( 0, nvar )
        kappa[2,] <- mu.central[2,]
        kappa[3,] <- mu.central[3,]
        for ( j in 1:nvar ) {
            n <- 3
            while ( n <= order.max ) {
                np1 <- n + 1
                nm1 <- n - 1
                total <- mu.raw[np1,j]
                for ( k in 1:nm1 ) {
                    km1 <- k - 1
                    total <- total -  choose(nm1,km1) * kappa[k+1,j] * mu.raw[n-k+1,j]
                }
                kappa[np1,j] <- total
                n <- n + 1
            }
        }
        return( kappa )
    }
    else if ( is.data.frame( mu.raw ) ) {
        order.max <- nrow( mu.raw ) - 1
        nvar <- ncol( mu.raw )
        kappa <- data.frame( matrix( nrow=order.max+1,ncol=nvar ) )
        kappa[1,] <- rep( 0, nvar )
        kappa[2,] <- mu.central[2,]
        kappa[3,] <- mu.central[3,]
        for ( j in 1:nvar ) {
            n <- 3
            while ( n <= order.max ) {
                np1 <- n + 1
                nm1 <- n - 1
                total <- mu.raw[np1,j]
                for ( k in 1:nm1 ) {
                    km1 <- k - 1
                    total <- total -  choose(nm1,km1) * kappa[k+1,j] * mu.raw[n-k+1,j]
                }
                kappa[np1,j] <- total
                n <- n + 1
            }
        }
        return( kappa )
    }
    else
        stop( "argument mu.raw is not a vector, matrix or data frame" )
    return( NULL )
}
