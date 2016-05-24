
library( RMTstat )

testdata <- new.env()
load( "testdata.rda", envir=testdata )

test.dtw <- function( beta ) {
    x <- get( "x", testdata )
    d <- get( paste( "dtw", beta, sep='' ), testdata )
    max( abs( dtw( x, beta ) - d ) ) < 5e-16
}

test.dtw.log <- function( beta ) {
    x  <- get( "x", testdata )
    d  <- dtw( x, beta )
    ld <- dtw( x, beta, log=TRUE )
    max( abs( log( d ) - ld ) ) < 5e-14
}

test.ptw <- function( beta ) {
    x <- get( "x", testdata )
    p <- get( paste( "ptw", beta, sep='' ), testdata )
    max( abs( ptw( x, beta ) - p ) ) < 5e-16
}

test.ptw.tail <- function( beta ) {
    x  <- get( "x", testdata )
    p  <- ptw( x, beta )
    pu <- ptw( x, beta, lower.tail=FALSE )
    
    identical( 1-p, pu )
}

test.ptw.log <- function( beta ) {
    x  <- get( "x", testdata )
    p  <- ptw( x, beta)
    lp <- ptw( x, beta, log=TRUE )
    max( abs( log( p ) - lp ) ) < 5e-14
}

test.ptw.log.tail <- function( beta ) {
    x  <- get( "x", testdata )
    p   <- ptw( x, beta)
    lpu <- ptw( x, beta, log=TRUE, lower.tail=FALSE )
    max( abs( log( 1-p ) - lpu ) ) < 5e-14
}

test.qtw <- function( beta ) {
    p <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.5, 
           0.80, 0.90, 0.95, 0.99, 0.995, 0.999)
    x <- qtw( p, beta )
    max( abs( ptw( x, beta ) - p ) ) < 5e-8
}

test.qtw.tail <- function( beta ) {
    p <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.5, 
           0.80, 0.90, 0.95, 0.99, 0.995, 0.999)
    x <- qtw( p, beta, lower.tail=FALSE )
    max( abs( ptw( x, beta, lower.tail=FALSE ) - p ) ) < 5e-8
}

test.qtw.log <- function( beta ) {
    p <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.5, 
           0.80, 0.90, 0.95, 0.99, 0.995, 0.999)
    x <- qtw( log( p ), beta, log.p=TRUE )
    max( abs( ptw( x, beta ) - p ) ) < 5e-8
}

test.qtw.log.tail <- function( beta ) {
    p <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.5, 
           0.80, 0.90, 0.95, 0.99, 0.995, 0.999)
    x <- qtw( log( p ), beta, lower.tail=FALSE, log.p=TRUE )
    max( abs( ptw( x, beta, lower.tail=FALSE ) - p ) ) < 5e-8
}

test.WishartSpikePar <- function() {
    n      <- 10
    p      <- 99
    ratio  <- p/n
    var    <- 0.28401
    spike  <- ( sqrt( ratio ) + 0.01 )*var
    center <- ( spike + var )*( 1 + ratio*( var/spike ) )
    scale  <- ( spike + var )*sqrt( 2*( 1 - ratio*( var/spike )^2 )/n )
    params <- WishartSpikePar( spike, n, p, var )
    
    identical( params, list( centering=center, scaling=scale ) )
}

tests <- list( function() { test.ptw( 1 ) }
             , function() { test.ptw( 2 ) }
             , function() { test.ptw( 4 ) }

             , function() { test.ptw.tail( 1 ) }
             , function() { test.ptw.tail( 2 ) }
             , function() { test.ptw.tail( 4 ) }
             
             , function() { test.ptw.log( 1 ) }
             , function() { test.ptw.log( 2 ) }
             , function() { test.ptw.log( 4 ) }

             , function() { test.ptw.log.tail( 1 ) }
             , function() { test.ptw.log.tail( 2 ) }
             , function() { test.ptw.log.tail( 4 ) }

             , function() { test.qtw( 1 ) }
             , function() { test.qtw( 2 ) }
             , function() { test.qtw( 4 ) }

             , function() { test.qtw.tail( 1 ) }
             , function() { test.qtw.tail( 2 ) }
             , function() { test.qtw.tail( 4 ) }

             , function() { test.qtw.log( 1 ) }
             , function() { test.qtw.log( 2 ) }
             , function() { test.qtw.log( 4 ) }
                          
             , function() { test.qtw.log.tail( 1 ) }
             , function() { test.qtw.log.tail( 2 ) }
             , function() { test.qtw.log.tail( 4 ) }

             , function() { test.dtw( 1 ) }
             , function() { test.dtw( 2 ) }
             , function() { test.dtw( 4 ) }
             
             , function() { test.dtw.log( 1 ) }
             , function() { test.dtw.log( 2 ) }
             , function() { test.dtw.log( 4 ) }

             , function() { test.WishartSpikePar() }
             )

test.all <- function() {
    n.tests <- length( tests )
    passed  <- rep( NA, n.tests)

    cat( "Running tests" )
    for( i in 1:n.tests ) {
        passed[ i ] <- tests[[i]]()
        cat( "." )
    }
    cat( "done.\n" )

    n.passed <- sum( passed )
    n.failed <- n.tests - n.passed

    if( n.failed > 0 ) {
        cat( "There were", n.failed, "failures:\n\n" )
        for( i in 1:n.tests ) {
            if( !passed[ i ] ) {
                print( tests[[ i ]] )
                cat( "\n" )
            }
        }
    }
    
    cat( "Tests run: ", n.tests, ", ",
         "Failures: ", n.tests-n.passed, "\n", sep='' )
    
    invisible( n.failed == 0 )
}

test.all()
