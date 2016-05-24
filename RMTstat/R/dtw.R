dtw <- function( x, beta=1, log=FALSE ) {
    if( beta != 1 && beta != 2 && beta != 4 )
        stop( "'beta' must be '1', '2', or '4'.")
        
    x <- as.vector( x )
    d <- rep( NA, length( x ) )
    
    below <- x < q0.xmin
    above <- x > q0.xmax
    idx   <- !below & !above    
    x     <- x[ idx ]

    if( !log ) {
        d[ below ] <- 0.0
        d[ above ] <- 0.0

        if( beta == 1 ) {
            d[ idx ] <- ( -q0sol( x, 4 ) + q0sol( x, 1 ) )*ptw( x, 1 )/2
        } else if( beta == 2 ) {
            d[ idx ] <- ptw( x, 2 ) * ( -q0sol( x, 4 ) )
        } else if( beta == 4 ) {
            d[ idx ] <- ( ( 1/( 2*sqrt( ptw( x, 2 ) ) ) ) 
                            * dtw( x, 2 ) 
                            * cosh( q0sol( x, 5 )/2 ) 
                          + sqrt( ptw( x, 2 ) ) 
                            * sinh( q0sol( x, 5 )/2 ) 
                            * (-q0sol( x, 1 )/2 ) )
        }
    } else {
        d[ below ] <- -Inf
        d[ above ] <- -Inf

        if( beta == 1 ) {
            d[ idx ] <- ( log( -q0sol( x, 4 ) + q0sol( x, 1 ) ) 
                          + ptw( x, 1, log.p=TRUE ) 
                          - log( 2 ) )
        } else if( beta == 2 ) {
            d[ idx ] <- ptw( x, 2, log.p=TRUE ) + log( -q0sol( x, 4 ) )
        } else if( beta == 4 ) {
            d[ idx ] <- log( dtw( x, 4, log=FALSE ) )
        }        
    }
    
    d   
}
