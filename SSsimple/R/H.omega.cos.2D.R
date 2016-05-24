H.omega.cos.2D <-
function(x, y, u.x, u.y, phs.x, phs.y) {
	n <- length(x)
	d <- length(u.x) * length(u.y)
	H <- matrix(NA, n, d)
	for(i in 1:n) {
		H[ i, ] <- as.vector( tcrossprod( cos(y[i]*u.y+phs.y) , cos(x[i]*u.x+phs.x) )  )	
	}
	return(H)
}

