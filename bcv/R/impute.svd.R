
impute.svd.check <- function( impute ) {
    function( x, k=min(n,p), tol=max(n,p)*1e-10, maxiter=100 ) {
        x  <- as.matrix( x )
        n  <- nrow( x )
        p  <- ncol( x )
        np <- min( n, p )
    
        if( is.complex( x ) )
            stop( "'x' cannot be complex." )
        if( !all( is.finite( x ) | is.na( x ) ) )
            stop( "'x' cannot contain any infinities or NaNs." )
        if( !( k >= 0 ) )
            stop( "'k' must be in the non-negative, not '", k, "'." )
        if( !( k <= np ) )
            stop( "'k' must be less than '", np, "' not '", k, "'." )
        if( !( tol > 0 ) )
            stop( "'tol' must be positive, not, '", tol, "'." )
        if( !( maxiter > 0 ) )
            stop( "'maxiter' must be a positive value, not '", 
                  maxiter, "'." )
    
        storage.mode( x ) <- "double"
    
        res <- impute( x, k, tol, maxiter )
        
        if( res$iter == maxiter )
            warning( "Did not converge in '", maxiter, "' iterations." )
            
        res
    }
}

impute.colMeans <- function( x ) {
    isna      <- is.na( x )
    x[ isna ] <- 0
    
    if( min( dim( x ) ) > 0 ) {
        s  <- drop( matrix( 1, 1, nrow( x ) ) %*% x )
        n  <- nrow( x ) - apply( isna, 2, sum )
        mu <- ifelse( n == 0, 0, s/n )
    } else {
        mu <- rep( 0, ncol( x ) )
    }
    
    mu
}

impute.svd.R.unchecked <- function( x, k, tol, maxiter )
{
    n <- nrow( x )
    p <- ncol( x )
    K <- seq_len( k )
      
    if (n == 0 || p == 0) {
        return( list( x=x, rss=0, iter=0 ) )
    }

    missing     <- is.na( x )
    missing.idx <- which( missing )
    missing.j   <- ( ( missing.idx-1 ) %/% n ) + 1

    mu                   <- impute.colMeans( x )
    missing.est          <- mu[ missing.j ]
    xhat0                <- x
    xhat0[ missing.idx ] <- missing.est
    
    rss0  <- Inf
    delta <- Inf
    iter  <- 0
    
    while( ( delta > tol ) && ( iter < maxiter ) ) {
        xhat0[ !missing ] <- x[ !missing ]
        xhat0.svd         <- svd( xhat0 )
        uhat              <- xhat0.svd$u[ ,K,drop=FALSE ] 
        dhat              <- xhat0.svd$d[ K ] 
        vhat              <- xhat0.svd$v[ ,K,drop=FALSE ]
        
        xhat1 <- uhat %*% ( dhat * t( vhat ) )
        
        rss1  <- sum( ( ( x-xhat1 )[ !missing ] )^2 )
        delta <- abs( rss0-rss1 )/( .Machine$double.eps + rss1 )
        
        iter  <- iter + 1
        xhat0 <- xhat1
        rss0  <- rss1
    }
        
    x[ missing ] <- xhat0[ missing ]
    
    list( x=x, rss=rss0, iter=iter )
}
impute.svd.R <- impute.svd.check( impute.svd.R.unchecked )

impute.svd.C.unchecked <- function( x, k, tol, maxiter )
{
    res <- .Call( "R_impute_svd", x, k, tol, maxiter )
    names( res ) <- c("x", "rss", "iter")
    res
}
impute.svd.C <- impute.svd.check( impute.svd.C.unchecked )

impute.svd <- impute.svd.C
