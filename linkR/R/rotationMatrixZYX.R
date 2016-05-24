rotationMatrixZYX <- function(t, t2=NULL, t3=NULL){

	if(length(t) == 1){
		t[2] <- t2
		t[3] <- t3
	}

	z <- matrix(c(cos(t[1]), -sin(t[1]), 0, sin(t[1]), cos(t[1]), 0, 0, 0, 1), 3, 3, byrow=TRUE)
	y <- matrix(c(cos(t[2]), 0, sin(t[2]), 0, 1, 0, -sin(t[2]), 0, cos(t[2])), 3, 3, byrow=TRUE)
	x <- matrix(c(1, 0, 0, 0, cos(t[3]), -sin(t[3]), 0, sin(t[3]), cos(t[3])), 3, 3, byrow=TRUE)

	z %*% y %*% x
}