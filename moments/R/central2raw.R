central2raw <- function( mu.central, eta )
{
    if ( is.vector( mu.central ) ) {
        if ( length( eta ) > 1 )
            stop( "argument eta has too many values" )
        np1 <- length( mu.central )
        n <- np1 - 1
        mu.raw <- rep( 0, np1 )
        mu.raw[1] <- 1
        j <- 2
        for ( k in 1:n ) {
            total <- 0
            for ( i in 0:k ) {
                total <- total + choose(k,i) * ( ( eta )^(k-i) ) * mu.central[i+1]
            }
            mu.raw[j] <- total
            j <- j + 1
        }
        return( mu.raw )
    }
    mu.raw <- NULL
    if ( is.matrix( mu.central ) ) {
        np1 <- nrow( mu.central )
        nrv <- ncol( mu.central )
        mu.raw <- matrix( nrow=np1, ncol=nrv )
    }
    if ( is.data.frame( mu.central ) ) {
        np1 <- nrow( mu.central )
        nrv <- ncol( mu.central )
        mu.raw <- data.frame( matrix( nrow=np1, ncol=nrv ) )
    }
    if ( is.null( mu.raw ) )
        stop( "argument mu.central is not a vector, matrix or data frame" )
    if ( length( eta ) != nrv )
        stop( "argument eta does not match the argument mu.central" )
    n <- np1 - 1
    for ( m in 1:nrv ) {
        mu.raw[1,m] <- 1
        eta.m <- eta[m]
        j <- 2
        for ( k in 1:n ) {
            total <- 0
            for ( i in 0:k ) {
                total <- total + choose(k,i) * ( ( eta.m )^(k-i) ) * mu.central[i+1,m]
            }
            mu.raw[j,m] <- total
            j <- j + 1
        }
    }
    return( mu.raw )
        
}
