ptw <- function( q, beta=1, lower.tail = TRUE, log.p = FALSE )
{
    if( beta != 1 && beta != 2 && beta != 4 )
        stop( "'beta' must be '1', '2', or '4'.")
        
    q <- as.vector( q )
    p <- rep( NA, length( q ) )
    
    below <- q < q0.xmin
    above <- q > q0.xmax
    idx   <- !below & !above    
    q     <- q[ idx ]

    if( !log.p ) {
        p[ below ] <- 0.0
        p[ above ] <- 1.0

        if( beta == 1 ) {
            p[ idx ] <- sqrt( ptw( q, 2 ) ) * exp( -q0sol( q, 5 )/2 )
        } else if( beta == 2 ) {
            p[ idx ] <- exp( -q0sol( q, 3 ) )
        } else if( beta == 4 ) {
            p[ idx ] <- sqrt( ptw( q, 2 ) ) * cosh( q0sol( q, 5 )/2 )
        }
        
        if( !lower.tail ) {
            p <- 1 - p
        }
        
    } else {
        if( lower.tail ) {
            p[ below ] <- -Inf
            p[ above ] <-  0

            if( beta == 1 ) {
                p[ idx ] <- 0.5 * ptw( q, 2, log.p=TRUE ) -q0sol( q, 5 )/2
            } else if( beta == 2 ) {
                p[ idx ] <- -q0sol( q, 3 )
            } else if( beta == 4 ) {
                p[ idx ] <- ( 0.5 * ptw( q, 2, log.p=TRUE )
                              + log( cosh( q0sol( q, 5 )/2 ) ) )
            }
        } else {
            p[ below ] <- 0
            p[ above ] <- -Inf
            p[ idx ]   <- log( ptw( q, beta, lower.tail=FALSE, log.p=FALSE ) )
        }
    }
    
    p
}
