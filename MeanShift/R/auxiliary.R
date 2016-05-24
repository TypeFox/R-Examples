gaussianKernel <- function( x ){
	
	## function to evaluate the asymmetric gaussian kernel	
	computeGaussianKernel <- function( y ){
	
		if( 0 <= y ){
		
			value <- 2 / 0.388 * dnorm( y / 0.388 )
		
		} else{
		
			value <- 0
		
		}
	
		return( value )
	
	}
	
	output <- sapply( x, computeGaussianKernel )
	
	return( output )
		
}


###

exponentialKernel <- function( x ){
	
	## function to evaluate the asymmetric exponential kernel	
	computeExponentialKernel <- function( y ){
	
		if( 0 <= y ){
		
			value <- dexp( y, rate=4.61 )
		
		} else{
		
			value <- 0
		
		}
	
		return( value )
	
	}
	
	output <- sapply( x, computeExponentialKernel )
	
	return( output )
		
}

###

cubicKernel <- function( x ){
	
	## function to evaluate the asymmetric cubic kernel	
	computeCubicKernel <- function( y ){
	
		if( 0 <= y && y<= 1 ){
		
			value <- 4 * ( 1 - y )^3
		
		} else{
		
			value <- 0
		
		}
	
		return( value )
	
	}
	
	output <- sapply( x, computeCubicKernel )
	
	return( output )
		
}

###

epanechnikovKernel <- function( x ){
	
	## function to evaluate the asymmetric Epanechnikov kernel	
	computeEpanechnikovKernel <- function( y ){
	
		if( 0 <= y && y<= 1 ){
		
			value <- 3 / 2 * ( 1 - y^2 )
		
		} else{
		
			value <- 0
		
		}
	
		return( value )
	
	}
	
	output <- sapply( x, computeEpanechnikovKernel )
	
	return( output )
		
}

###

distanceFunction <- function( x, y ){
	
	## function to compute the standard euclidean distance
	output <- sqrt( sum( ( x - y )^2 ) )
	
	return( output )
	
}

###

connectedComponents <- function( X, tol.epsilon=1e-3 ){

	N <- ncol( X )
	
	## initialize components matrix
	C <- X
	
	## initialize components vector
	labels <- vector( mode="integer", length=N )
	
	K <- 1 
	labels[1] <- 1
	C[,1] <- X[,1]
	
	# pb <- txtProgressBar( min=0, max=N, style=3 )
	
	## efficient connected component algorithm
	for( n in 2:N ){
		
		assigned <- FALSE
				
		for( k in 1:K ){
			
			distance <- distanceFunction( X[,n], C[,k] )
			
			if( distance < tol.epsilon ){
				
				labels[n] <- k
				assigned <- TRUE
				break
				
			}
			
		}
		
		if( !assigned ){
			
			K <- K + 1
			labels[n] <- K
			C[,K] <- X[,n]
			
		}
		
		# setTxtProgressBar( pb, n )
		
	}
	
	C <- as.matrix( C[,1:K] )
	colnames( C ) <- paste( "mode", 1:K, sep="" )
	
	labels <- as.integer( labels )
	
	output <- list( components=C, labels=labels )
	
	# close( pb )
	
	message( "\nThe algorithm found ", as.character( K ),
	" clusters.\n")
	
	return( output )
		
}