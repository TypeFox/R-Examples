print.predict.mex <- function( x, ... ){
    if ( is.R() ) stdev <- function( x ) sqrt( var( x ) )

    print( x$call, ... )

    cv <- names( x$data$simulated )[ 1 ]
    dn <- paste( "E(", dimnames(x$data$simulated)[[ 2 ]] ,"|", cv , ">Q",100*x$pqu,")", sep="" )

    if (!is.null(x$replicates)){
        cat( "\nResults from", length( x$replicates ), "bootstrap runs.\n" )
        res <- t( sapply( x$replicates , function ( x ) apply( x, 2, mean ) ) )
        m <- apply( res, 2, mean )
        s <- apply( res, 2, stdev )
        res <- matrix( c( m, s ), byrow=TRUE, nrow=2 )
        dimnames( res ) <- list( c( "mean", "se" ), dn )
    }
    else {
        res <- matrix(apply(x$data$simulated, 2, mean), nrow=1)
        dimnames(res) <- list("mean", dn)
    }

    cat( paste( "\nConditioned on ", cv, " being above its ", 100*x$pqu, "th percentile.\n\n", sep = "" ) )
    print(res, ...)
}