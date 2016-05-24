dltInverse <- function(cal.coeff, coor.3d){
	# Translated from Matlab function dlt_inverse written by Ty Hedrick

	# IF INPUT IS VECTOR CONVERT TO SINGLE ROW MATRIX
	if(!is.matrix(coor.3d)){coor.3d <- matrix(coor.3d, 1, 3)}

	# SOLVE FOR 2D COORDINATES GIVEN CALIBRATION COEFFICIENTS AND 3D COORDINATES
	m = matrix(NA, nrow(coor.3d), 2)
	m[, 1] <- (coor.3d[, 1]*cal.coeff[1] + coor.3d[, 2]*cal.coeff[2] + coor.3d[, 3]*cal.coeff[3] + cal.coeff[4]) / (coor.3d[, 1]*cal.coeff[9] + coor.3d[, 2]*cal.coeff[10] + coor.3d[, 3]*cal.coeff[11] + 1)
	m[, 2] <- (coor.3d[, 1]*cal.coeff[5] + coor.3d[, 2]*cal.coeff[6] + coor.3d[, 3]*cal.coeff[7] + cal.coeff[8]) / (coor.3d[, 1]*cal.coeff[9] + coor.3d[, 2]*cal.coeff[10] + coor.3d[, 3]*cal.coeff[11] + 1)

	m
}
