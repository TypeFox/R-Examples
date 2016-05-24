

round.fold <- function( n, k ) {
    if (n <= 1)
        stop( "n should be greater than 1" )
    
    k     <- round( k )
    kvals <- unique( round( n/seq_len( floor( n/2 ) ) ) )
    temp  <- abs( kvals - k )
    if (!any( temp == 0 ) ) 
        k <- kvals[ temp == min( temp ) ][1]
    k
}

choose.sets <- function( n, k ) {
    if (k < 1)
        stop( "k should be positive" )
    if (k > n)
        stop( "k should be less than or equal to n" )
    
    f   <- ceiling( n/k )
    s   <- sample( rep( seq_len( k ), f ), n )
    n.s <- table( s )

    if ( length( n.s ) != k )
        choose.sets(n,k)
    else
        s
}
