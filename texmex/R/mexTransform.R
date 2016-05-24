`mexTransform` <-
function( x ,
		  method = "mixture",
		  divisor = "n+1",
		  na.rm = TRUE,
      margins = "laplace"){

	if ( !is.element( method, c( "mixture", "empirical" ) ) )
		stop( "method should be either 'mixture' or 'empirical'" )

	if ( divisor == "n" )
		divisor <- dim( x$data )[[ 1 ]]
	else if ( divisor == "n+1" )
		divisor <- dim( x$data )[[ 1 ]] + 1
	else stop( "divisor can be 'n' or 'n+1'" )

	transFun <- function( i , x , mod , th, divisor, method ){
		x <- x[ , i ]
		mod <- mod[[ i ]]
		th <- th[ i ]

		ox <- order( x )
		names( x ) <- 1:length( x )
	
		x <- sort( x )
		run <- rle( x )
		p <- cumsum( run$lengths )/ divisor
		p <- rep( p, run$lengths )
	
		p <- p[ order( as.character( names( x ) ) ) ]
		x <- x[ order( as.character( names( x ) ) ) ]
	
	
		Femp <- p
		if ( method == "mixture" ){
			sigma <- exp( mod$coefficients[ 1 ] )
			xi <- mod$coefficients[ 2 ]
			
			Para <- ( 1 + xi * ( x - th) / sigma ) ^ ( -1 / xi )
			Para <- 1 - mean( x > th ) * Para
			res <- ifelse( x <= th , Femp , Para )
		}
		else res <- Femp

		res[ ox ] <- sort( res )
		res
	} # Close transfun
	
	res <- sapply( 1:( dim( x$data )[[ 2 ]] ), transFun,
				      x = x$data, mod = x$models, th = x$mth,
				      divisor = divisor, method=method
				     )
	dimnames( res ) <- list( NULL, names( x$models ) )

    if (margins == "gumbel"){
    	x$transformed <- -log( -log( res ) )
    }
    else if (margins == "laplace"){
		#x$transformed <- sign(res - .5) * log(1 - 2*abs(res - .5))

		x$transformed <- ifelse(res < .5,  log(2 * res), -log(2 * (1 - res) ))
    }
    else { stop("margins need to be gumbel or laplace") }
	invisible(x)
}

