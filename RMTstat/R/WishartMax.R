
WishartMaxPar <- (function() {
    mu <- function( n,p ) {
        n.sqrt <- sqrt( n )
        p.sqrt <- sqrt( p )
        res    <- ( n.sqrt + p.sqrt )^2
        res
    }
    
    sigma <- function( n,p ) {
        n.sqrt <- sqrt( n )
        p.sqrt <- sqrt( p )
        res    <- ( n.sqrt + p.sqrt )*( 1/n.sqrt + 1/p.sqrt )^( 1/3 )
        res
    }

    mu.real <- function( n,p ) {
        mu( n-1/2,p-1/2 )
    }
    
    sigma.real <- function( n,p ) {
        sigma( n-1/2,p-1/2 )
    }
    
    alpha <- function( n,p ) {
        1/( 1 + ( mu( n-1/2,p+1/2 )/mu( n+1/2,p-1/2 ) )
                * sqrt( sigma( n-1/2,p+1/2 )/sigma( n+1/2,p-1/2 ) ) )
    }
    
    mu.cplx <- function( n,p ) {
        a   <- alpha( n,p )
        res <- mu( n-1/2,p+1/2 )*a + mu( n+1/2,p-1/2 )*( 1-a )
        res
    }
    
    sigma.cplx <- function( n,p ) {
        a   <- alpha( n,p )
        res <- sigma( n-1/2,p+1/2 )*a + sigma( n+1/2,p-1/2 )*( 1-a )
        res
    }
    
    function( ndf, pdim, var=1, beta=1 ) {
        n <- ndf
        p <- pdim
        
        if( beta == 1 ) {
            m <- mu.real( n,p )
            s <- sigma.real( n,p )
        } else if( beta == 2 ) {
            m <- mu.cplx( n,p )
            s <- sigma.cplx( n,p )
        } else {
            stop( "`beta' must be 1 or 2, not `", beta, "'")
        }
            
        center <- var*( m/n )
        scale  <- var*( s/n )
    
        list( centering=center, scaling=scale )
    }
})()

dWishartMax <- function( x, ndf, pdim, var=1, beta=1, log = FALSE ) {
    params <- WishartMaxPar( ndf, pdim, var, beta )
    x.tw   <- ( x - params$centering )/( params$scaling )
    d.tw   <- dtw( x.tw, beta, log )
    d <- if (log)
        d.tw - log( params$scaling )
    else d.tw / ( params$scaling )
    d
}

pWishartMax <- function( q, ndf, pdim, var=1, beta=1, 
                      lower.tail = TRUE, log.p = FALSE ) {
    params <- WishartMaxPar( ndf, pdim, var, beta )
    q.tw   <- ( q - params$centering )/( params$scaling )
    p      <- ptw( q.tw, beta, lower.tail, log.p )
    p
}

qWishartMax <- function( p, ndf, pdim, var=1, beta=1,
                      lower.tail = TRUE, log.p = FALSE ) {
    params <- WishartMaxPar( ndf, pdim, var, beta )
    q.tw   <- qtw( p, beta, lower.tail, log.p )
    q      <- params$centering + q.tw*( params$scaling )
    q
}

rWishartMax <- function( n, ndf, pdim, var=1, beta=1 ) {
    params <- WishartMaxPar( ndf, pdim, var, beta )
    x.tw   <- rtw( n, beta )
    x      <- params$centering + x.tw*( params$scaling )
    x
}
