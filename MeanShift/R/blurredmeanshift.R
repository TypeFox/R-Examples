blurringMeanShiftOperator <- function( X, h=1, kernel="epanechnikovKernel" ){
	
	n.curves <- ncol( X )
	
	## compute distances
	distances <- as.matrix( dist( t( X ), diag=TRUE, upper=TRUE ) )
	
	## scale by bandwidth
	scaled.distances <- distances / h
	
	## evaluate kernel
	kernel <- get( kernel )
	kernel.values <- matrix( kernel( scaled.distances ), nrow=n.curves,
	ncol=n.curves ) 
	
	## weights denominators
	total.sum <- colSums( kernel.values )
	
	## weights
	kernel.weights <- kernel.values / total.sum

	## update
	new.X <- X%*%t( kernel.weights )
	
	output <- new.X
	
	return( new.X )
	
}

blurringMeanShiftAlgorithm <- function( X, h=NULL,
kernel="epanechnikovKernel", tol.stop=1e-6, max.iter=100 ){
	
	if( is.null( h ) ){
		
		h <- quantile( dist( t( X ) ), 0.3 )
		
	}
	
	close.enough <- FALSE
	
	old.X <- X
	
	iter <- 0
	not.converged <- FALSE
	
	## while the largest update corresponds to a shift
	## larger than 'tol.stop' and while number of iterations
	## is smaller than 'max.iter'
	while( !close.enough ){
		
		## apply blurring mean-shift operator and update
		iter <- iter + 1
		
		new.X <- blurringMeanShiftOperator( X=old.X, h=h, kernel=kernel )
		
		distance <- max( sqrt( colSums( old.X - new.X )^2 ) )
		
		old.X <- new.X
		
		close.enough <- ( distance < tol.stop )
		
		if( iter >= max.iter ){
			
			not.converged <- TRUE
			break
			
		}
		
	}
	
	if( not.converged ){
		
		if( kernel == "epanechnikovKernel"){
			
			warning( "Reached maximum number of iterations (", 
			as.character( max.iter),"). The algorithm ",
			"didn't converge. Try increasing max.iter." )
			
		} else{

			warning( "Reached maximum number of iterations (", 
			as.character( max.iter),"). The algorithm ",
			"didn't converge. Try kernel=\"epanechnikovKernel\"." )
			
		}
		
	} else{

		message( "Blurring mean-shift algorithm ran successfully.\n")
			
	}
	
	return( new.X )
	
}

bmsClustering <- function( X, h=NULL, kernel="epanechnikovKernel",
tol.stop=1e-6, max.iter=100, tol.epsilon=1e-3 ){
	
	# minimal input checking
	X <- as.matrix( X )
	max.iter <- as.integer( max.iter )
	
	if( ncol( X ) <= 1 ){
		
		message( "The input matrix X has only one column: ",
		"returning input.")
		return( X )
	}

	if( !is.element( kernel, paste( c( "epanechnikov", "cubic", 
	"gaussian", "exponential"), "Kernel", sep="" ) ) ){
		
		stop( "Invalid kernel name.")
		
	}
	
	if( !is.null( h ) && h <= 0 ){
		
		stop( "The bandwidth must be strictly positive." )
				
	}
	
	if( max.iter <= 0 ){
		
		stop( "The maximum number of iterations must be a positive ",
		"integer." )
		
	}
	
	if( tol.stop <= 0 || tol.epsilon <= 0 ){
		
		stop( "All tolerances must be strictly positive.")
		
	}
	
	## run blurring mean-shift algorithm
	message( "\nRunning blurring mean-shift algorithm...\n" )
	
	blurring.mean.shift.algorithm <- blurringMeanShiftAlgorithm( X=X,
	h=h, kernel=kernel, tol.stop=tol.stop, max.iter=max.iter )
	
	## find connected components
	message( "Finding clusters..." )
	output <- connectedComponents( X=blurring.mean.shift.algorithm,
	tol.epsilon=tol.epsilon )
	
	invisible( output )

}