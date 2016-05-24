contr.bayes.ordered <- function ( n, scores = 1:n, contrasts = TRUE )
{
    make.poly <- function( n, scores ) {
        y   <- scores - mean( scores )
        X   <- outer( y, seq_len( n ) - 1, "^" )
        QR  <- qr( X )
        z   <- QR$qr
        z   <- z *( row( z ) == col( z ) )
        raw <- qr.qy( QR, z )
        Z   <- sweep( raw, 2, apply( raw, 2, function( x ) sqrt( sum( x^2 ) ) ), "/" )
        colnames( Z ) <- paste( "^", 1:n - 1, sep="" )
        Z
    }
    if ( is.numeric( n ) && length( n ) == 1 ) { levs <- 1:n }
    else {
        levs <- n
        n <- length( levs )
    }
    if ( n < 2 ) {
        stop( gettextf( "contrasts not defined for %d degrees of freedom", n - 1 ), domain = NA ) 
    }
    if ( n > 95 ) {
        stop( gettextf( "orthogonal polynomials cannot be represented accurately enough for %d degrees of freedom", n-1 ), domain = NA ) 
    }
    if ( length( scores ) != n ) {
        stop( "'scores' argument is of the wrong length" )
    }
    if ( !is.numeric( scores ) || any( duplicated( scores ) ) ) {
        stop("'scores' must all be different numbers")
    }
    contr <- make.poly( n, scores )
    if ( contrasts ) {
        dn <- colnames( contr )
        dn[2:min( 4, n )] <- c( ".L", ".Q", ".C" )[1:min( 3, n-1 )]
        colnames( contr ) <- dn
        contr[, , drop = FALSE]
    }
    else {
        contr[, 1] <- 1
        contr
    }
}

contr.bayes.unordered <- function(n, base = 1, contrasts = TRUE) {
    if( is.numeric( n ) && length( n ) == 1) {
        if( n > 1 ) { levs <- 1:n }
        else stop( "not enough degrees of freedom to define contrasts" )
    } 
    else {
        levs <- n
        n <- length( n )
    }
    contr <- array( 0, c(n, n), list( levs, levs ) )
    diag( contr ) <- 1
    if( contrasts ) {
        if( n < 2 ) { stop( gettextf( "contrasts not defined for %d degrees of freedom", n - 1 ), domain = NA ) }
        if( base < 1 | base > n ){ stop( "baseline group number out of range" ) }
        contr <- contr[, , drop = FALSE]
    }
    contr
}
