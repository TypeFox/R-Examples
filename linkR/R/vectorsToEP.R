vectorsToEP <- function(u, v, h='right'){

	axis <- cprod(u, v, h)
	angle <- avectors(u, v)
	R <- tMatrixEP(axis, angle)
	euler_angles <- rotationMatrixToEP(R)

	euler_angles
}