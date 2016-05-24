threshold_slope <- function( pfit, xfit, thresh = 0.5 ) {
#
# The function finds the approximate value of x (=x_th) for which the value of the 
# estimated psychometric function is equal to 'thresh' and the approximate value of 
# slope in x_th.
# 
# INPUT
# 
# pfit - estimated values of the psychometric function
# xfit - stimulus levels at which the function was estimated
#
# OPTIONAL INPUT
# 
# thresh - criterion level at which to estimate threshold; default is 0.5
#
# OUTPUT
#
# Object with 2 elements:
# 	x_th  - estimated threshold
# 	slope - estimated value of slope, i.e. derivative of pfit at x_th

# MAIN PROGRAM
# First two arguments are mandatory
    if( missing("pfit") || missing("xfit") ) {
        stop("Check input. First 2 arguments are mandatory");
    }
    
    if( !is.double( thresh ) || length( thresh ) > 1 ) {
        stop( "Threshold level must be scalar" );
    }
    
    if( ( thresh < 0 ) || (thresh > 1 )) {
        stop( "Threshold level must be betwen 0 and 1" );
    }
   
# Check that the input variables match
    
if ( length( pfit ) != length( xfit ) )
    stop( 'Length of fitted values pfit must be the same as length of xfit' );
end


    value <- NULL;
# threshold
    value$x_th <- xfit[ which( abs( pfit - thresh ) == min( abs( pfit - thresh ) ) ) ];

    if( length( value$x_th ) > 1 ) {
# if there are many point for the same threshold value, then function is flat
# in this point and slope=0
        value$slope <- 0;
        value$x_th <- mean( value$x_th );
    }
    else {
# slope
        ind <- which( xfit == value$x_th );
        value$slope <- ( pfit[pmin( ind + 1, length( pfit ) )] - 
                       pfit[pmax( ind - 1, 1 )] ) /
                     ( xfit[pmin( ind + 1, length( xfit ) )] -
                       xfit[pmax( ind - 1, 1 )] );
    }
    return( value );
}