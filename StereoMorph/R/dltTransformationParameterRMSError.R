dltTransformationParameterRMSError <- function(p, coor.2d, nx, ny, sx, sy=NULL, p.fixed=NULL){

	#t_iter <<- t_iter + 1

	# GENERATE 3D CALIBRATION POINTS USING TRANSFORMATION PARAMETERS
	coor.3d <- transformPlanarCalibrationCoordinates(tpar=c(p.fixed, p), nx=nx, ny=ny, sx=sx, sy=sy)

	# FIND COEFFICIENTS AND DLT INVERSE ERROR
	dlt_coefficients <- dltCoefficients(coor.3d=coor.3d, coor.2d=coor.2d)
	#cat('\tMean RMS Error: ', round(mean(dlt_coefficients$rmse), 6), "\n", sep="")

	# SET MINIMIZATION ERROR
	rmse_error <- mean(dlt_coefficients$rmse, na.rm=TRUE)

	rmse_error
}