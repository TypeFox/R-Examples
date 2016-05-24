scaleX <- function ( x, a=min( x, na.rm=TRUE ), b=max( x, na.rm=TRUE ), u, v )
{
###
### This function returns the values in the vector x that have been mapped
### from the interval [a,b] to the interval [u,v]
###
### Parameters
### x = a numerical vector of values to be mapped
### na.rm = a boolean value that determines what to do with missing values
### a = a numerical value for the lower bound of the domain interval (default is min(x) )
### b = a numerical value for the upper bound of the domain interval (default is max(x) )
### u = a numerical value for the lower bound of the target interval
### v = a numerical value for the upper bound of the target interval
###
    ###
    ### parameter validation section
    ###
    if ( missing( x ) )
        stop( "x is missing" )
    if ( !is.numeric( x ) )
        stop( "x does not contain numerical values" )
    if ( !is.vector( x ) )
        stop( "x is not a vector" )
    if ( !is.numeric( a ) || !is.vector( a ) )
        stop( "domain lower bound parameter a must be a numeric value" )
    if ( length( a ) > 1 )
        stop( "domain lower bound parameter a is not a scalar value" )
    if ( !is.numeric( b ) || !is.vector( b ) )
        stop( "domain upper bound parameter b must be a numeric value" )
    if ( length( b ) > 1 )
        stop( "domain upper bound parameter b is not a scalar value" )
    if ( a > b )
        stop( "domain lower bound a is greater than upper bound b" )
    if ( missing( u ) )
        stop( "target lower bound parameter u is missing" )
    if ( missing( v ) )
        stop ( "target upper bound parameter v is missing" )
    if ( !is.numeric( u ) || !is.vector( u ) )
        stop( "lower bound u must be a numeric value" )
    if ( length( u ) > 1 )
        stop( "target lower bound parameter u is not a scalar value" )
    if ( !is.numeric( v ) || !is.vector( v ) )
        stop( "upper bound v must be a numeric value" )
    if ( length( v ) > 1 )
        stop( "target upper bound parameter v is not a scalar value" )
    if ( u > v )
        stop( "target lower bound u is greater than upper bound v" )
    ###
    ### default result
    ###
    result <- NULL    
    ###
    ### the target [-Inf,Inf] is unchanged
    ###
    if ( ( u == -Inf ) && ( v == Inf ) )
        result <- x
    ###
    ### the semi-infinite target interval [-Inf,v] is a shift
    ###
    if ( ( u == -Inf ) && is.finite( v ) )
        result <- x + v - b
    ###
    ### the semi-infinite target interval [u,Inf] is a shift
    ###
    if ( is.finite( u ) || v == Inf )
        result <- x + u - a
    ###
    ### the finite target interval requires a shift and scaling
    ###
    if ( is.finite( u ) && is.finite( v ) )
        result <- u + ((v - u) * (x - a) / (b - a))
    ###
    ### attach the interval attributes to the result
    ###
    attr( result, "a" ) <- a
    attr( result, "b" ) <- b
    attr( result, "u" ) <- u
    attr( result, "v" ) <- v
    return( result )
}    
