qqtrunc <- function( x, spec, a=-Inf, b=Inf, title="Truncated Distribution Q-Q Plot", 
    xlabel="Theoretical Quantiles", ylabel="Sample Quantiles", ... )
{
###
### This function produces a q-q plot of sample quantiles against theoretical
### quantiles for a truncated random variable
###
### Arguments
### x   = a numeric vector of sample values
### min = a numeric value for the lower bound of the uniform distribution
### max = a numeric value for the upper bound of the uniform distribution
###
    plot( qtrunc(ppoints(x), spec, a=a, b=b, ...), 
        sort(x), main=title,
        xlab=xlabel, ylab=ylabel ) 
    abline( 0,1 )
    return( invisible( 0 ) )
}
