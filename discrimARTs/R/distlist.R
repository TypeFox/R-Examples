## Each distlist function must return 2 valid distribution (pdf) functions, 
## each taking a vector of quantiles and 2 parameters.
## Distribution functions taking >2 parameters can be used, 
## but can only be optimized over 2 (positional) parameters;
## The remaining params should be encapsulated within the function calls.
## mix.distlist.facing.gamma provides a good example for users to extend

mix.distlist.norm <- function() {
    ## construct a list containing 2 normal pdf functions
    ## for mk.mle's distlist arg
    distlist <- list( dist1 = dnorm, dist2 = dnorm )
    return(distlist)
}

mix.distlist.facing.gamma <- function(lower, upper) {
    ## construct a list containing 2 facing gamma pdf functions
    ## requires lower and upper bound of data
    ## for mk.mle's distlist arg
    ## NOTE: shape and scale are passed by position: 
    ## shape is par1, scale is par2
    distlist <- list()
    ## forward
    distlist$dist1 <- function(x, .shape, .scale) {
        dgamma(x-lower, shape=.shape, scale=.scale)
    }
    ## reverse
    distlist$dist2 <- function(x, .shape, .scale) {
        dgamma(upper-x, shape=.shape, scale=.scale)
    }
    return(distlist)
}
