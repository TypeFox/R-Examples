createErrorPlots <- function(cal.coeff, corners, nx, sq.size.num, sq.size.units = '', 
	file = ''){

	# TEST CALIBRATION ACCURACY AGAINST OPTIM SET
	dlt_test <- dltTestCalibration(cal.coeff, corners, nx, paste0(sq.size.num, ' ', sq.size.units))
	
	if(is.null(file)) return(dlt_test)

	# ADD DIRECTORY SYMBOL IF NOT EMPTY
	if(file != '') if(!grepl('/$', file)) file <- paste0(file, '/')

	# CREATE PLOTS
	pdf(file=paste0(file, 'Epipolar error by aspect.pdf'))
	par(mar=c(6.1, 4.1, 3.1, 2.1))
	boxplot(dlt_test$epipolar.error, ylab='Epipolar error (pixels)', xlab='', las=2)
	mtext('Aspect', side=1, line=4)
	dev.off()

	pdf(file=paste0(file, 'Point-to-point length error by aspect.pdf'))
	par(mar=c(6.1, 4.1, 3.1, 2.1))
	boxplot(dlt_test$ipd.error, ylab=paste0('Point-to-point length error (', sq.size.units, ')'), xlab='', las=2)
	abline(h=0, lty=3, col=gray(0.5))
	mtext('Aspect', side=1, line=4)
	dev.off()

	pdf(file=paste0(file, 'Point-to-point length percent error by aspect.pdf'))
	par(mar=c(6.1, 4.1, 3.1, 2.1))
	boxplot((dlt_test$ipd.error / dlt_test$pair.dist)*100, ylab=paste0('Point-to-point length error (%)'), xlab='', las=2)
	abline(h=0, lty=3, col=gray(0.5))
	mtext('Aspect', side=1, line=4)
	dev.off()

	#pdf(file=paste0(file, 'Adjacent point-to-point length error by aspect.pdf'))
	#par(mar=c(6.1, 4.1, 3.1, 2.1))
	#boxplot(dlt_test$adj.pair.ipd.error, ylab=paste0('Adjacent point-to-point length error (', sq.size.units, ')'), xlab='', las=2)
	#mtext('Aspect', side=1, line=4)
	#dev.off()
	
	pdf(file=paste0(file, 'Point-to-point length, error v length.pdf'))
	plot(c(dlt_test$pair.dist), c(dlt_test$ipd.error), 
		xlab=paste0('Point-to-point length (', sq.size.units, ')'), 
		ylab=paste0('Point-to-point length error (', sq.size.units, ')'), type='n')
	abline(h=0, lty=3, col=gray(0.5))
	#abline(a=0, b=0.05, lty=2, col=gray(0.5))
	#abline(a=0, b=-0.05, lty=2, col=gray(0.5))
	points(c(dlt_test$pair.dist), c(dlt_test$ipd.error), 
		xlab=paste0('Point-to-point length (', sq.size.units, ')'), 
		ylab=paste0('Point-to-point length error (', sq.size.units, ')'))
	dev.off()

	pdf(file=paste0(file, 'Epipolar error.pdf'))
	hist <- hist(c(dlt_test$epipolar.error), main=paste0('Epipolar error\n(', dim(corners)[3], ' aspects)'),
		xlab='Epipolar error (pixels)')
	ee_mean <- mean(dlt_test$epipolar.error, na.rm=TRUE)
	ee_sd <- sd(dlt_test$epipolar.error, na.rm=TRUE)
	ee_max <- max(dlt_test$epipolar.error, na.rm=TRUE)
	abline(v=ee_mean, col='purple', lty=2, lwd=2)
	abline(v=ee_sd, col='green', lty=3, lwd=2)
	abline(v=ee_max, col='orange', lty=4, lwd=2)
	legend(x=mean(hist$breaks), y=max(hist$counts), 
		legend=c(
			paste0('Mean: ', round(ee_mean, 2), ' px'),
			paste0('SD: ', round(ee_sd, 2), ' px'),
			paste0('Max: ', round(ee_max, 2), ' px')
		), 
		lty=2:4, lwd=2, col=c('purple', 'green', 'orange'),
		xjust=0.5, yjust=1, bty='n')
	dev.off()

	pdf(file=paste0(file, 'Point-to-point length error.pdf'))
	hist <- hist(c(dlt_test$ipd.error), main=paste0('Point-to-point length error\n(', dim(corners)[3], ' aspects)'),
		xlab=paste0('Point-to-point length error (', sq.size.units, ')'))
	ipd_mean <- mean(dlt_test$ipd.error, na.rm=TRUE)
	ipd_abs_mean <- mean(abs(dlt_test$ipd.error), na.rm=TRUE)
	ipd_sd <- sd(abs(dlt_test$ipd.error), na.rm=TRUE)
	ipd_max <- max(abs(dlt_test$ipd.error), na.rm=TRUE)
	abline(v=ipd_mean, col='blue', lty=2, lwd=2)
	abline(v=ipd_abs_mean, col='purple', lty=3, lwd=2)
	abline(v=ipd_sd, col='green', lty=4, lwd=2)
	abline(v=ipd_max, col='orange', lty=5, lwd=2)
	legend(x=min(hist$breaks), y=max(hist$counts), 
		legend=c(
			paste0('Mean: ', round(ipd_mean, 3), ' ', sq.size.units),
			paste0('Abs mean: ', round(ipd_abs_mean, 3), ' ', sq.size.units),
			paste0('Abs SD: ', round(ipd_sd, 3), ' ', sq.size.units),
			paste0('Abs Max: ', round(ipd_max, 3), ' ', sq.size.units)
		), 
		lty=2:5, lwd=2, col=c('blue', 'purple', 'green', 'orange'),
		xjust=0, yjust=1, bty='n')
	dev.off()

	pdf(file=paste0(file, 'Point-to-point length percent error.pdf'))
	p_error <- (c(dlt_test$ipd.error) / c(dlt_test$pair.dist))*100
	hist <- hist(p_error, main=paste0('Point-to-point length percent error\n(', dim(corners)[3], ' aspects)'),
		xlab=paste0('Point-to-point length percent error (%)'))
	ipd_mean <- mean(p_error, na.rm=TRUE)
	ipd_abs_mean <- mean(abs(p_error), na.rm=TRUE)
	ipd_sd <- sd(abs(p_error), na.rm=TRUE)
	ipd_max <- max(abs(p_error), na.rm=TRUE)
	abline(v=ipd_abs_mean, col='blue', lty=2, lwd=2)
	abline(v=ipd_mean, col='purple', lty=3, lwd=2)
	abline(v=ipd_sd, col='green', lty=4, lwd=2)
	abline(v=ipd_max, col='orange', lty=5, lwd=2)
	legend(x=min(hist$breaks), y=max(hist$counts), 
		legend=c(
			paste0('Mean: ', round(ipd_mean, 3), ' %'),
			paste0('Abs mean: ', round(ipd_abs_mean, 3), ' %'),
			paste0('Abs SD: ', round(ipd_sd, 3), ' %'),
			paste0('Abs Max: ', round(ipd_max, 3), ' %')
		), 
		lty=2:5, lwd=2, col=c('blue', 'purple', 'green', 'orange'),
		xjust=0, yjust=1, bty='n')
	dev.off()
	

	### Reconstruct points and align along major axes
	# Remove NA corners
	corners <- corners[, , rowSums(is.na(corners[1,1,,])) == 0, ]
	
	# Make 3D corner array
	corners_3d <- array(NA, dim=c(dim(corners)[1], 3, dim(corners)[3]), 
		dimnames=list(NULL, c('x', 'y', 'z'), dimnames(corners)[[3]]))
	rec_errors <- matrix(NA, nrow=dim(corners)[1], ncol=dim(corners)[3],
		dimnames=list(NULL, dimnames(corners)[[3]]))

	# FILL 3D CORNER ARRAY
	for(aspect in 1:dim(corners)[3]){

		# RECONSTRUCT CORNERS INTO 3D
		dlt_recon <- dltReconstruct(cal.coeff, corners[, , aspect, ])
		rec_errors[, aspect] <- dlt_recon$rmse
		corners_3d[, , aspect] <- dlt_recon$coor.3d
	}

	# Find centers of each checkerboard
	corners_3d_centers <- t(apply(corners_3d, 3, 'colMeans'))

	# Find the corner centroid
	corners_centroid <- colMeans(corners_3d_centers)
	
	# Center all corners about centroid
	corners_3d <- corners_3d - array(matrix(corners_centroid, nrow=dim(corners_3d)[1], ncol=3, byrow=TRUE), dim=dim(corners_3d))

	# Find centers of translated checkerboards
	corners_3d_centers <- t(apply(corners_3d, 3, 'colMeans'))

	# Find major axis of 3D points
	pca <- princomp(corners_3d_centers)
	prin_comp <- rbind(cprod_SM(pca$loadings[, 'Comp.1'], pca$loadings[, 'Comp.2']), pca$loadings[, 'Comp.2'], pca$loadings[, 'Comp.1'])

	# Align the principal components with z,y,x
	RM <- tMatrixDC_SM(prin_comp, diag(3))
	
	# Apply rotation matrix
	for(aspect in 1:dim(corners_3d)[3]) corners_3d[, , aspect] <- corners_3d[, , aspect] %*% RM
	
	# Plot reconstruction error as a function of position along the major axes
	pdf(file=paste0(file, 'Reconstruction error v position along major axes.pdf'), width=5.5, height=7.5)
	par(mfrow=c(3,1), mar=c(4,4,1,1))
	plot(corners_3d[, 1, ], rec_errors, ylab='Reconstruction error (px)', xlab=paste0('Position along first major axis (', sq.size.units, ')'))
	plot(corners_3d[, 2, ], rec_errors, ylab='Reconstruction error (px)', xlab=paste0('Position along second major axis (', sq.size.units, ')'))
	plot(corners_3d[, 3, ], rec_errors, ylab='Reconstruction error (px)', xlab=paste0('Position along third major axis (', sq.size.units, ')'))
	dev.off()

	# Find interpoint distance errors
	ipd_pos <- array(NA, dim=c(floor(dim(corners_3d)[1]/2), 3, dim(corners_3d)[3]), dimnames=list(NULL, c('x','y','z'), dimnames(corners_3d)[[3]]))
	ipd_error <- matrix(NA, nrow=floor(dim(corners_3d)[1]/2), ncol=dim(corners_3d)[3], dimnames=list(NULL, dimnames(corners_3d)[[3]]))

	for(aspect in 1:dim(corners_3d)[3]){
		ipd_list <- findInterpointDistanceError(coor.3d=corners_3d[, , aspect], nx=nx, ny=dim(corners_3d)[1]/nx, sq.size=sq.size.num)
		ipd_pos[, , aspect] <- ipd_list$ipd.pos
		ipd_error[, aspect] <- ipd_list$ipd.error
	}

	# Plot interpoint distance error as a function of position along the major axes
	pdf(file=paste0(file, 'Point-to-point length error v position along major axes.pdf'), width=5.5, height=7.5)
	par(mfrow=c(3,1), mar=c(4.5,4.5,1,1))
	ylab <- paste0('Point-to-point length error (', sq.size.units, ')')
	plot(ipd_pos[, 1, ], ipd_error, ylab=ylab, xlab=paste0('Position along first major axis (', sq.size.units, ')'))
	abline(h=0, lty=2)
	plot(ipd_pos[, 2, ], ipd_error, ylab=ylab, xlab=paste0('Position along second major axis (', sq.size.units, ')'))
	abline(h=0, lty=2)
	plot(ipd_pos[, 3, ], ipd_error, ylab=ylab, xlab=paste0('Position along third major axis (', sq.size.units, ')'))
	abline(h=0, lty=2)
	dev.off()
	
	dlt_test
}