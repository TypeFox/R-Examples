dltCoefficientRMSError <- function(p, coor.2d){

	#c_iter <<- c_iter + 1

	# GENERATE 3D CALIBRATION POINTS USING TRANSFORMATION PARAMETERS
	dlt_reconstruct <- dltReconstruct(matrix(p, nrow=11, ncol=dim(coor.2d)[3]), coor.2d)
	#cat('Mean RMS Error: ', mean(dlt_reconstruct$rmse), '\n')

	#if(print.progress && c_iter %% 10 == 0) cat(c_iter, ' ')

	mean(dlt_reconstruct$rmse)
}