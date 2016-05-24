
cv.svd.wold.check <- function( cv.svd ) {
    function( x, k=5, maxrank=20, tol=1e-4, maxiter=20 ) {
        x <- as.matrix( x )
        n <- nrow( x )
        p <- ncol( x )
    
        if (is.complex( x ))
            stop ("x cannot be complex")
        if (n*p < 2)
            stop ("x should have at least two elements")
        if (k > n*p || k <= 1)
            stop("k outside allowable range")
        if (maxrank < 0 || maxrank > min(n,p))
            stop("maxrank outside allowable range")
    
        storage.mode( x ) <- "double"
    
        k.o <- k; k <- round.fold( n*p, k );
        if (k != k.o) 
            warning("k has been set to ", k)
    
        sets  <- choose.sets( n*p, k )
        
        msep <- cv.svd( x, k, maxrank, tol, maxiter, sets )
        colnames( msep ) <- 0:maxrank
        
        res   <- list( call=match.call(), k=k, maxrank=maxrank, 
                       msep=msep, sets=sets )
        class( res ) <- c("cvsvd_wold", "cvsvd")
        res
    }
}

cv.svd.wold.R.unchecked <- function( x, k, maxrank, tol, maxiter, sets ) {
    msep <- matrix( NA, k, maxrank+1 )
    
    for( j in seq_len( k ) ) {
        train  <- sets != j
        test   <- !train
        
        xtrain <- x
        xtrain[ test ] <- NA
        
        for( rank in seq( 0, maxrank ) ) {
            xhat <- suppressWarnings(
                        impute.svd( xtrain, rank, tol, maxiter )$x )
                        
            err  <- mean( ( xhat[ test ] - x[ test ] )^2 )
            msep[ j, rank+1 ] <- err
        }
    }
    
    msep
}
cv.svd.wold.R <- cv.svd.wold.check( cv.svd.wold.R.unchecked )

cv.svd.wold.C.unchecked <- function( x, k, maxrank, tol, maxiter, sets ) {
    msept <- .Call( "R_cv_svd_wold", x, k, maxrank, tol, maxiter, 
                     as.integer( sets-1 ) )    
    msep  <- t( msept )
    msep
}
cv.svd.wold.C <- cv.svd.wold.check( cv.svd.wold.C.unchecked )

cv.svd.wold <- cv.svd.wold.C
