cprod_SM <- function(u, v, h = 'right'){

	# CROSS PRODUCT OF TWO VECTORS IN THREE-DIMENSIONAL SPACE
	r <- c(u[2]*v[3] - u[3]*v[2], u[3]*v[1] - u[1]*v[3], u[1]*v[2] - u[2]*v[1])

	r
}