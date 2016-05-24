


#***************************************************************
# RRUM parametrization
.rrum.param <- function( delta.summary , q.matrix ){
	#---
	#  RRUM parametrization
	#  log( P(X=1) ) = b0 + b1*alpha1 + b2 * alpha2 
	#  RRUM:
	#  P(X=1) = pi * r1^( 1- alpha1) * r2^(1-alpha2)
	#  => log( P(X=1) ) = log[ pi * r1 * r2 * r1^(-alpha1) * r2^(-alpha2) ]
	#                   = log( pi ) + log(r1) + log(r2) + -log(r1)*alpha1 + -log(r2) * alpha2
	#  => b1 = -log(r1) and r1 = exp( -b1 )
	#  => log(pi) = b0 + b1 + b2 and pi = exp( b0 + b1 + b2 )
	I <- nrow(q.matrix)
	K <- ncol(q.matrix)
	rrum.params <- matrix( NA , I , K+1 )
	rownames(rrum.params) <- delta.summary[ delta.summary$partype == 0 , "item" ]
	colnames(rrum.params) <- c( "pi" , paste( "r_", colnames(q.matrix) , sep="") )
	for (ii in 1:I){
		# ii <- 2
		d.ii <- delta.summary[ delta.summary$itemno == ii , ]
		rrum.params[ii,"pi"] <- exp( sum( d.ii$est ) )
		rrum.params[ ii , which( q.matrix[ii,]==1) +1 ] <- exp( - d.ii$est[-1] )
				}
	return( rrum.params )
        }