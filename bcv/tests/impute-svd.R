
require( bcv )

impute.svd.R <- get( "impute.svd.R", 
                     env=parent.env( environment( impute.svd ) ) )


check <- function( actual, test ) {
    failed   <- FALSE
    expected <- suppressWarnings( 
                    impute.svd.R( test$x, test$k, test$tol, test$maxiter ) )
    
    if( !all( abs( expected$x - actual$x ) 
              < 
              1e-12 + 1e-8 * abs( expected$x ) ) ) {

        failed <- TRUE
    }
    else if( !( abs( expected$rss - actual$rss ) 
                <
                1e-12 + 1e-8 * abs( expected$rss ) ) ) {
                    
        failed <- TRUE            
    }
    else if (!( expected$iter == actual$iter ) ){

        failed <- TRUE
    } 

    if (failed) {
        cat( "Expected:\n")
        print( expected )
        
        cat( "Actual:\n")
        print( actual )
    }
    
    failed
}

rdim <- function( size ) { 
    floor( runif( 1, 0, size ) ) 
}

rmatrix <- function( m,n ) { 
    matrix( rnorm( m*n ), m, n )
}

rproportion <- function() { 
    sample( c(0, 1, runif( 1 )), 1, prob=c(0.1, 0.1, 0.8) ) 
}

rmissing <- function( m, n, prop ) {
    mn       <- m * n;
    nmissing <- round( prop*mn )
    
    if (mn > 0) {
        idx <- sample( seq_len( mn ), nmissing )
    } else {
        idx <- c()
    }
    
    idx
}

rmatrix.partial <- function( m, n ) {
    x      <- rmatrix( m, n )
    p      <- rproportion()
    idx    <- rmissing( m, n, p )
    x[idx] <- NA
    x
}

rrank <- function( m, n ) {
    sample( 0:min( m, n ), 1)
}

rtol <- function() {
    sample( c(1e-5, 1e-6, 1e-7, 1e-8), 1 )
}

rmaxiter <- function() {
    sample( c(1, 10, 50, 100), 1 )
}

rtest <- function( size ) {
    m       <- rdim( size )
    n       <- rdim( size )
    x       <- rmatrix.partial( m, n )
    k       <- rrank( m, n )
    tol     <- rtol()
    maxiter <- rmaxiter()

    list( x=x, k=k, tol=tol, maxiter=maxiter )
}

sizes <- function( ntests, each=2 ) {
    if( ntests > 0 )
        rep(4 + sqrt( seq( 0, ntests/each, length=ceiling( ntests/each ) ) )
                    , each=each )[ 1:ntests ]
    else
        c()
}


ntests <- 1000
s <- sizes( ntests )
set.seed( 0 )

nsuc <- 0
for (size in s)
{
    cat( '.' )
    test <- rtest( size )
    
    actual <- suppressWarnings( 
                  impute.svd( test$x, test$k, test$tol, test$maxiter ) )
    if( !check( actual, test ) )
        nsuc <- nsuc + 1
}

if( nsuc == ntests ) {
    cat("Passed", nsuc, "tests!\n")
} else {
    cat("Not all tests passed.\n")
}

