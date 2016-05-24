qqunif <- function( x, u.min=0, u.max=1, title="Uniform Q-Q Plot", xlabel="Theoretical Quantiles",
    ylabel="Sample Quantiles" )
{
###
### This function produces a q-q plot of sample quantiles against theoretical
### quantiles for a uniformly distributed random variable
###
### Arguments
### x = a numeric vector of sample values
### u.min = a numeric value for the lower bound of the uniform distribution
### u.max = a numeric value for the upper bound of the uniform distribution
###
    return( invisible( plot( qunif(ppoints(x), min=u.min, max=u.max), sort(x), main=title,
        xlab=xlabel, ylab=ylabel ) ) )
}
