ghermite.h.quadrature <- function( functn, rule, mu=0, lower=-Inf, upper=Inf, 
weighted=TRUE,  ... )
{
###
### This function evaluates the integral of the function functn
### between lower and upper using the weight and abscissa values specified
### in the rule data frame.  The rule corresponds to an order n
### Hermite polynomial, weight function and interval [-Inf,Inf]
### Lower bound is finite and upper bound is finite.
###
### Parameters
### functn   = an R function which should take a numeric argument x and
###            possibly some parameters.  The function returns a
###            numerical vector value for the given argument x.
### rule     = a data frame containing the order n quadrature rule
### mu       = numeric value for the polynomial parameter
### lower    = a scalar lower bound of the integral
### upper    = a scalar lower found of the integral.
### weighted = a boolean value which if true includes the weight function in the integrand
### ...      = other arguments passed to the function functn.
###
    if ( !is.function( functn ) )
        stop( "functn argument is not an R function" )
    if ( !is.data.frame( rule ) )
        stop( "rule argument is not a data frame" )
    if ( weighted ) {
        ff <- 
            if ( length( list( ... ) ) && length( formals( functn ) ) > 1 )
                function( x ) functn( x, ... )
            else
                functn
    }
    else {
        ff <- 
            if ( length( list( ... ) ) && length( formals( functn ) ) > 1 )
                function( x ) { functn( x, ... ) / ghermite.h.weight( x, mu ) }
            else
                function( x ) { functn( x ) / ghermite.h.weight( x, mu ) }
 }
    if ( !is.infinite( lower ) )
        stop( "lower limit l is not infinite" )
    if ( !is.infinite( upper ) )
        stop( "upper limit u is not infinite" )
    if ( lower == upper )
        return( 0 )
        s <- if ( ( upper == Inf ) && ( lower == -Inf ) )
                +1
            else
                -1
    x <- rule$x
    w <- rule$w
    return( s * sum( w * ff( s * x ) ) )
}
