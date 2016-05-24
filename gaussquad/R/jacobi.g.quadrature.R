jacobi.g.quadrature <- function( functn, rule,  p=1, q=1, lower=0, upper=1, 
weighted=TRUE, ... )
{
###
### This function evaluates the integral of the function functn
### between lower and upper using the weight and abscissa values specified
### in the rule data frame.  The rule corresponds to an order n
### Jacobi polynomial, weight function and interval [0,1]
### Lower bound is finite and upper bound is finite.
###
### Parameters
### functn   = an R function which should take a numeric argument x and
###            possibly some parameters.  The function returns a
###            numerical vector value for the given argument x.
### rule     = a data frame containing the order n quadrature rule
### p        = first polynomial parameter
### q        = second polynomial parameter
### lower    = a scalar lower bound of the integral
### upper    = a scalar lower found of the integral.
### weighted = a boolean value which if true includes the weight function in the integrand
### ...      = other arguments passed to the function functn.
###
    if ( !is.function( functn ) )
        stop( "functn argument is not an R function" )
    if ( is.infinite( lower ) )
        stop( "lower bound is infinite" )
    if ( is.infinite( upper ) )
        stop( "lower bound is infinite" )
    if ( weighted ) {
        ff <- 
            if ( length( list( ... ) ) && length( formals( functn ) ) > 1 )
                function( x, p, q ) { functn( x, ... ) }
            else
                function( x, p, q ) { functn( x ) }
	}
    else {
        ff <- 
            if ( length( list( ... ) ) && length( formals( functn ) ) > 1 )
                function( x, p, q ) { functn( x, ... ) / jacobi.g.weight( x, p, q ) }
            else
                function( x, p, q ) { functn( x ) / jacobi.g.weight( x, p, q ) }
    }
    lambda <- upper - lower
    mu <-     lower
    y <- lambda * rule$x + mu
    w <- rule$w
    return( lambda * sum( w * ff(y) ) )
}
