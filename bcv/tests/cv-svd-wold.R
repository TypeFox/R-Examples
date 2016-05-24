
require( bcv )

cv.svd.wold.R <- get( "cv.svd.wold.R.unchecked", 
                      env=parent.env( environment( cv.svd.wold ) ) )


check <- function( actual, test ) {
    failed   <- FALSE
    expected <- suppressWarnings( 
                    cv.svd.wold.R( test$x, length( unique( actual$sets ) ), 
                                   test$maxrank, test$tol, 
                                   test$maxiter, actual$sets ) )
    
    if (!identical( dim( actual$msep ), dim( expected ) ) ) {
        failed <- TRUE
    } else if( !all( abs( expected - actual$msep )
                     < 
                     1e-12 + 1e-8 * abs( expected ) ) ) {

        failed <- TRUE
    }

    if (failed) {
        cat( "Expected:\n")
        print( expected )
        
        cat( "Actual:\n")
        print( actual )

        browser()
    }
    
    failed
}

rdim <- function( size ) { 
    floor( runif( 1, 2, size ) ) 
}

rmatrix <- function( m,n ) { 
    matrix( rnorm( m*n ), m, n )
}

rfold <- function( m, n ) {
    sample( seq( 2, m*n ), 1, )
}

rmaxrank <- function( m, n ) {
    sample( seq_len( min( m, n ) ), 1 )
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
    x       <- rmatrix( m, n )
    k       <- rfold( m, n )
    maxrank <- rmaxrank( m, n )
    tol     <- rtol()
    maxiter <- rmaxiter()

    list( x=x, k=k, maxrank=maxrank, tol=tol, maxiter=maxiter )
}

sizes <- function( ntests, each=2 ) {
    if( ntests > 0 )
        rep(4 + sqrt( seq( 0, ntests/each, length=ceiling( ntests/each ) ) )
                    , each=each )[ 1:ntests ]
    else
        c()
}

main <- function () {
    ntests <- 100
    s <- sizes( ntests )
    set.seed( 0 )

    nsuc <- 0
    for (size in s)
    {
        cat( '.' )
        test <- rtest( size )
    
        actual <- suppressWarnings( 
                      cv.svd.wold( test$x, test$k, test$maxrank, test$tol, 
                                   test$maxiter ) )
        if( !check( actual, test ) )
            nsuc <- nsuc + 1
    }

    if( nsuc == ntests ) {
        cat("Passed", nsuc, "tests!\n")
    } else {
        cat("Not all tests passed.\n")
    }
}

main()

