raw2central <- function( mu.raw )
{
    if ( is.matrix( mu.raw ) )
        return( apply( mu.raw, 2, raw2central ) )
    else if ( is.vector( mu.raw ) ) {
        np1 <- length( mu.raw )
        n <- np1 - 1
        mu.central <- rep( 0, np1 )
        mu.central[1] <- 1
        eta <- mu.raw[2]
        j <- 2
        for ( k in 1:n ) {
            total <- 0
            for ( i in 0:k ) {
                total <- total + choose(k,i) * ( (- eta )^(k-i) ) * mu.raw[i+1]
            }
            mu.central[j] <- total
            j <- j + 1
        }
        return( mu.central )
    }
    else if ( is.data.frame( mu.raw ) )
        return( sapply( mu.raw, raw2central ) )
    else
        return( raw2central( as.vector( mu.raw ) ) )
}
