
MarchenkoPasturPar <- function( ndf=NA, pdim=NA, var=1, svr=ndf/pdim ) {
    gamma          <- svr
    
    inv.gamma.sqrt <- sqrt( 1/gamma )
    a              <- var*( 1 - inv.gamma.sqrt )^2
    b              <- var*( 1 + inv.gamma.sqrt )^2
    
    list( lower=a, upper=b )
}

dmp <- function( x, ndf=NA, pdim=NA, var=1, svr=ndf/pdim, log = FALSE ) {
    gamma  <- svr

    params <- MarchenkoPasturPar( ndf, pdim, var, svr )
    a      <- params$lower
    b      <- params$upper
    
    if( !log ) {
        # we have to handle +/- zero carefully when gamma=1
        d <- ifelse( gamma == 1 & x == 0 & 1/x > 0, Inf, 
                 ifelse( x <= a | x >= b, 0,
                     suppressWarnings(
                         gamma/( 2*pi*var*x )
                         *
                         sqrt( ( x-a )*( b-x ) ) 
                         ) ) )
    } else {
        d <- ifelse( gamma == 1 & x == 0 & 1/x > 0, Inf,
                 ifelse( x <= a | x >= b, -Inf,
                     suppressWarnings(
                         log( gamma ) 
                         - 
                         ( log( 2 ) + log( pi ) + log( var ) + log( x ) )
                         +
                         0.5*log( x-a ) 
                         + 
                         0.5*log( b-x )
                         ) ) )
    }
    d
}

pmp <- function( q, ndf=NA, pdim=NA, var=1, svr=ndf/pdim, 
                 lower.tail = TRUE, log.p = FALSE ) {
    params <- MarchenkoPasturPar( ndf, pdim, var, svr )
    a <- params$lower
    b <- params$upper
    
    f <- function( x ) {
             dmp( x, ndf, pdim, var, svr )
         }
    
    if( lower.tail ) {
        p <- ifelse( q <= a, 0, 
                 ifelse( q >= b, 1,
                     integrate( f, a, q )$value 
                 )
             )
        p <- ifelse( svr < 1 && q >= 0, p + (1 - svr), p )
    } else {
        p <- ifelse( q >= b, 0,
                 ifelse( q <= a, min( 1, svr ),
                     integrate( f, q, b )$value
                 )
             )
        p <- ifelse( svr < 1 && q <= 0, p + (1 - svr), p )
    }
    
    res <- ifelse( log.p, log( p ), p )
    res
}
pmp <- Vectorize( pmp )

qmp <- function( p, ndf=NA, pdim=NA, var=1, svr=ndf/pdim, 
                 lower.tail = TRUE, log.p = FALSE ) {
    p      <- ifelse( log.p, exp( p ), p )
    p      <- ifelse( lower.tail, p, 1-p )
    params <- MarchenkoPasturPar( ndf, pdim, var, svr )

    q      <- NULL
    if ( p <= 0 ) {
        q <- ifelse( svr <= 1, -0, parmams$lower )
    } else if ( p >= 1 ) {
        q <- params$upper
    }

    if( svr < 1 ) {
        if( p < 1 - svr ) {
            q <- -0
        } else if( p == 1 - svr ) {
            q <- 0
        }
    }

    if( is.null( q ) ) {
        F <- function( x ) {
                 pmp( x, ndf, pdim, var, svr ) - p
             }
        q <- uniroot( F, interval=c(params$lower, params$upper) )$root
    }
    
    q
}
qmp <- Vectorize( qmp )

rmp <- function( n, ndf=NA, pdim=NA, var=1, svr=ndf/pdim ) {
    u <- runif( n )
    qmp( u, ndf, pdim, var, svr )
}
