
nr.numdiff <- function( ll0 , ll1 , ll2 , h , eps = 10^(-10) ){ 
    d1 <- ( ll1 - ll2  ) / ( 2 * h )    # negative sign?
	# second order derivative
	# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
	d2 <- ( ll1 + ll2 - 2*ll0 ) / h^2	
 	d2 <- ifelse( abs(d2) < eps , eps , d2 )			
	delta.change <- - d1 / d2
    return( delta.change )
        }