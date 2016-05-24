rotateBody <- function(m, p, v, a){

	if(is.vector(m)) m <- matrix(m, nrow=1, ncol=length(m))

	# IF INPUT MATRIX HAS NO ROWS, RETURN EMPTY MATRIX
	if(dim(m)[1] == 0) return(m)

	# GET TRANSLATION MATRIX
	tm <- matrix(p, nrow=nrow(m), ncol=ncol(m), byrow=TRUE)

	# MAKE ROTATION POINT ORIGIN
	m <- m - tm

	# FIND ROTATION MATRIX
	RM <- tMatrixEP(v=v, a=a)
	
	# ROTATE POINTS
	m <- m %*% RM

	# MOVE POINTS BACK TO ORIGINAL POSITION
	m <- m + tm
	
	m
}