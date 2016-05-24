orthogonalProjectionToLine <- function(p, l1 = NULL, l2 = NULL){

	# FINDS POSITION OF POINT X ON A LINE AT THE MINIMUM DISTANCE FROM POINT X TO THE INPUT POINT P
	# Adapted from : http://paulbourke.net/geometry/pointline/

	# GET LINE PARAMETERS IF L1 IS LIST
	if(is.list(l1)){
		# GET LINE POINTS IF IN l1 LIST
		if('l1' %in% names(l1)){l2 <- l1$l2;l1 <- l1$l1}
	
		# GET LINE POINTS IF m, b
		if('m' %in% names(l1)){l2 <- c(1, l1$m + l1$b);l1 <- c(0, l1$b)}
	
		# GET LINE POINTS IF a, b, c
		if('a' %in% names(l1)){l2 <- c(1, -(l1$a + l1$c)/l1$b);l1 <- c(0, -l1$c/l1$b)}
	}

	# TREAT POINT AS MATRIX FOR FUNCTION BUT IF INPUT POINT IS VECTOR RETURN POINT AS VECTOR
	if(is.vector(p)){p <- matrix(p, nrow=1, ncol=length(p));return_vector <- TRUE}else{return_vector <- FALSE}

	# IF INPUT IS 2D POINT(S) AND l1, l2, ADD ZERO THIRD DIMENSION
	return_2d <- FALSE
	if(ncol(p) == 2){
		p <- cbind(p, rep(0, nrow(p)))
		l1 <- c(l1, 0)
		l2 <- c(l2, 0)
		return_2d <- TRUE		
	}

	r <- matrix(NA, nrow=nrow(p), ncol=ncol(p))

	for(i in 1:nrow(p)){
		# IF L2 ARGUMENT IS NULL THEN TREAT L1 AS VECTOR AND GET SECOND POINT ON LINE THROUGH ORIGIN
		if(is.null(l2)) l2 <- 2*l1
	
		# CHECK THAT POINTS DEFINING LINE ARE NOT COINCIDENT
		if(sum((l2 - l1)^2) == 0) return(NA)
	
		# SOLVE FOR POSITION OF POINT ON LINE
		u <- sum((p[i, ] - l1)*(l2 - l1)) / sum((l2 - l1)^2)
		r[i, ] <- l1 + u*(l2 - l1)
	}

	# RETURN AS 2D VECTOR
	if(return_vector && return_2d) return(as.vector(r)[1:2])

	# RETURN AS 3D VECTOR
	if(return_vector && !return_2d) return(as.vector(r))

	# RETURN AS 2D MATRIX
	if(!return_vector && return_2d) return(r[, 1:2])

	# RETURN AS 3D MATRIX
	r
}