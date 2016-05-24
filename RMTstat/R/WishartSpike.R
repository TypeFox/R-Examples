
WishartSpikePar <- function( spike, ndf=NA, pdim=NA, var=1, beta=1 ) {
    ratio  <- pdim/ndf
    above  <- spike > sqrt( ratio )*var 
    center <- ifelse( !above, NA,
                  ( spike + var )*( 1 + ratio*( var/spike ) ) )
    scale  <- ifelse( !above, NA,
                  ( ( spike + var )
                    * sqrt( (2/beta)*( 1 - ratio*( var/spike )^2 ) ) 
                    / sqrt( ndf ) ) )
    
    list( centering=center, scaling=scale )
}

WishartSpikeVecPar <- function( spike, ndf=NA, pdim=NA, var=1, beta=1 ) {
    ratio <- pdim/ndf
    above <- spike > sqrt( ratio )*var 
    cor2  <- ifelse( !above, NA,
                ( ( spike^2 - ratio*var^4 )
                  / ( spike*( spike + ratio*var ) ) 
                  / sqrt( ndf ) ) )
    
    cor2
}

dWishartSpike <- function( x, spike, ndf=NA, pdim=NA, var=1, beta=1,
                            log = FALSE ) {
    params <- WishartSpikePar( spike, ndf, pdim, var, beta )
    d      <- dnorm( x, mean=params$centering, sd=params$scaling, log=log )
    d
}

pWishartSpike <- function( q, spike, ndf=NA, pdim=NA, var=1, beta=1,
                            lower.tail = TRUE, log.p = FALSE ) {
    params <- WishartSpikePar( spike, ndf, pdim, var, beta )
    p      <- pnorm( q, mean=params$centering, sd=params$scaling, lower.tail, log.p )
    p
}                                

qWishartSpike <- function( p, spike, ndf=NA, pdim=NA, var=1, beta=1,
                            lower.tail = TRUE, log.p = FALSE ) {
    params <- WishartSpikePar( spike, ndf, pdim, var, beta )
    q      <- qnorm( q, mean=params$centering, sd=params$scaling, lower.tail, log.p )
    q
}                                

rWishartSpike <- function( n, spike, ndf=NA, pdim=NA, var=1, beta=1 ) {
    params <- WishartSpikePar( spike, ndf, pdim, var, beta )
    x      <- qnorm( n, mean=params$centering, sd=params$scaling )
    x
}
