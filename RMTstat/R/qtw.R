qtw <- function( p, beta=1, lower.tail = TRUE, log.p = FALSE ) {
    q <- uniroot( function( x ) { ptw( x, beta, lower.tail, log.p ) - p }, 
             interval=c(-10.000000000000001, 6.000000000000001 ), 
             tol=1e-8 )$root
    q
}

qtw <- Vectorize( qtw )
