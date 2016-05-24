pointNormalOnLine <- function(pt, l1, l2=NULL){

	# FINDS POSITION OF POINT X ON A LINE AT THE MINIMUM DISTANCE FROM POINT X TO THE INPUT POINT
	# Adapted from : http://paulbourke.net/geometry/pointline/

	return_vector <- FALSE
	if(is.vector(pt)){
		pt <- matrix(pt, nrow=1, ncol=length(pt))
		return_vector <- TRUE
	}

	r <- matrix(NA, nrow=nrow(pt), ncol=ncol(pt))

	for(i in 1:nrow(pt)){

		# IF L2 ARGUMENT IS NULL THEN TREAT L1 AS VECTOR AND GET SECOND POINT ON LINE THROUGH ORIGIN
		if(is.null(l2)) l2 <- c(0,0,0)
	
		# CHECK THAT POINTS DEFINING LINE ARE NOT COINCIDENT
		if(sum((l2 - l1)^2) == 0) return(NA)
	
		# SOLVE FOR POSITION OF POINT ON LINE
		u <- sum((pt[i, ] - l1)*(l2 - l1)) / sum((l2 - l1)^2)
		r[i, ] <- l1 + u*(l2 - l1)
	}

	if(return_vector) return(as.vector(r))
	r
}