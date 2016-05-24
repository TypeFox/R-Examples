mix.synthetic.normal <- function(N=200, mix.prob=0.5, mu1, sd1, mu2, sd2) {
    ## Function to draw random deviates from mixture of normals
    ## Check parameters
    if ( any( c(mu1, sd1, mu2, sd2) <0 ) ) stop('Distribution parameters must be positive')
    if ( mix.prob <0 || mix.prob>1)  stop('mix.prob must be between 0 and 1')
    ## Is sample from lower or upper?
    .samples <- runif(N) < mix.prob
    ## number of samples from upper distribution 
    N.lower <- sum( .samples)
    ## Sample from lower distributions
    .samples[.samples==TRUE] <- rnorm( N.lower, mu1, sd1)
    ## Sample the rest from lower distributions
    .samples[.samples==FALSE] <- rnorm( N-N.lower, mu2, sd2)
    if (any(.samples<0) ) stop( 'Negative values in mixture.  Only positive values are supported for discrimARTs')
    return(.samples)
}

mix.synthetic.facing.gamma <- function(N=200, mix.prob=0.5, lower, upper, shape1, scale1, shape2, scale2) {
    ## draw random deviates from mixture of facing gammas
    ## Only keeping samples between lower and upper bounds
    ## Check parameters
    if ( any( c(lower, upper, shape1, scale1, shape2, scale2) <0 ) ) stop('Distribution parameters must be positive')
    if ( mix.prob <0 || mix.prob>1)  stop('mix.prob must be between 0 and 1')
    ## Is sample from lower or upper?
    .samples <- runif(N) < mix.prob
    ## number of samples from upper distribution 
    N.lower <- sum( .samples)
    ## Sample from lower distributions
    .samples[.samples==TRUE] <- rgamma( N.lower, shape=shape1, scale=scale1) + lower
    ## Sample the rest from lower distributions
    .samples[.samples==FALSE] <- upper - rgamma( N-N.lower, shape=shape2, scale=scale2)
    ## only keep samples within bounds
    .inbounds <- (.samples >= lower | .samples >= upper)
    .samples <- .samples[.inbounds]
    #if (any(.samples<0) ) stop( 'Negative values in mixture.  Only positive values are supported for discrimARTs')
    return(.samples)
}
