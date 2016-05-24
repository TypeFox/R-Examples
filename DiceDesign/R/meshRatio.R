meshRatio <- function(design){
	#----------------------------------
	# Compute the meshratio criterion
	# For a regular mesh, ratio=1
	# input : design of n experiments
	# Example : Meshratio(matrix(runif(40),20,2))
	#----------------------------------

	X <- as.matrix(design)
	n <- dim(X)[1]
	dimension <- dim(X)[2]

	if ( n < dimension ){
    		stop('Warning : the number of points is lower than the dimension')
	}

	# To check the experimental region
	if ( min(X)<0 || max(X)>1 ){
		warning("The design is rescaling into the unit cube [0,1]^d.")
		M <- apply(X,2,max)
		m <- apply(X,2,min)
		for (j in 1:dim(X)[2]){
			X[,j] <- (X[,j]-m[j])/(M[j]-m[j])
		}	
	}
	
   	DistanceMax  <- -1.0E30
   	DistanceMin  <-  1.0E30
  
   	for (i in 1:(n-1)) {
      	DistMin  <- 1.0E30
      	DistMax  <- -1.0E30
      	for (k in 1 : n){
           		if (i != k){   # if the differs from current point
				Dist <- 0
				for (j in 1 : dimension){
					Dist  <- Dist + (X[i,j] -X[k,j])*(X[i,j] - X[k,j])
               		}
				if (Dist > DistMax){ DistMax <- Dist; }
					if (Dist < DistMin){ DistMin <- Dist; }
          			}
      		}
     			if (DistanceMax < DistMin){
				DistanceMax <- DistMin
			}
     			if (DistanceMin > DistMin){
			DistanceMin <- DistMin
			}  
   		}
	ratio <- sqrt(DistanceMax/DistanceMin)
   	return(ratio)
}