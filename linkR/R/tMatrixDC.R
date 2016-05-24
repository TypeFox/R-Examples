tMatrixDC <- function(a, b=NULL){

	# Direction cosines method
	# a is global coordinate system unit vectors, b is local coordinate system unit vectors
	# The transformation matrix returned will transform points in the b coordinate system to the a coordinate system
	
	# IF B IS NULL USE IDENTITY MATRIX AS LOCAL COORDINATE SYSTEM
	if(is.null(b)) b <- diag(nrow(a))

	r <- matrix(NA, 3, 3)
	r[1, ] <- c(sum(a[1,]*b[1,]), sum(a[2,]*b[1,]), sum(a[3,]*b[1,]))
	r[2, ] <- c(sum(a[1,]*b[2,]), sum(a[2,]*b[2,]), sum(a[3,]*b[2,]))
	r[3, ] <- c(sum(a[1,]*b[3,]), sum(a[2,]*b[3,]), sum(a[3,]*b[3,]))

	return(r)
}