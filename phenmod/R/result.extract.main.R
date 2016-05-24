result.extract.main <- function(mask.grid, result.grid, 
			model="pim", interpolate=TRUE, 
			silent=FALSE, withOutliers=FALSE){
	outliers <- result.grid$outlier.bb + result.grid$outlier.lc
	outliers.na <- which(is.na(outliers)==TRUE)
	outliers[outliers.na] <- rep(0, length(outliers.na))

	if (!silent){ cat("Extracting doy's calculated by model:\n") }
	if (model == "pim"){
	values.model <- result.extract.sub(mask.grid, result.grid$doy.bb.pim, 
				result.grid$gk4.x, result.grid$gk4.y, outliers=outliers,
				silent=silent, withOutliers=withOutliers)$values
	} else {
		if (model == "tsm"){
			values.model <- result.extract.sub(mask.grid, result.grid$doy.bb.tsm, 
				result.grid$gk4.x, result.grid$gk4.y, outliers=outliers,
				silent=silent, withOutliers=withOutliers)$values
		} else {
			cat("Invalid model name!\n")	
			return(NULL)
		}
	}
	if (!silent){ cat("Extracting observed DoY's:\n") }
	values.observed <- result.extract.sub(mask.grid, result.grid$doy.bb.observed, 
				result.grid$gk4.x, result.grid$gk4.y, outliers=outliers,
				silent=silent, withOutliers=withOutliers)$values
	
	
	# calculate difference of calculated doy's to observed doy's
	values.dif <- values.model - values.observed

	if (!interpolate) {
		if (!silent){ cat("Mask values.. ",sep="") }
		# mask values of doy.model
		values.model <- result.extract.mask(mask.grid, values.model)

		# mask values of doy.observed
		values.observed <- result.extract.mask(mask.grid, values.observed)

		# mask values of doy.dif.tsm
		values.dif <- result.extract.mask(mask.grid, values.dif)
		if (!silent){ cat("Done!\n") }
	} else {
		if (!silent){ cat("Interpolate values:\n",sep="") }
		# interpolate values
		values.model <- result.extract.interpolate(mask.grid=mask.grid, 
					values=values.model, alt=mask.grid$alt, 
					x=mask.grid$x, y=mask.grid$y)

		values.observed <- result.extract.interpolate(mask.grid=mask.grid, 
					values=values.observed, alt=mask.grid$alt, 
					x=mask.grid$x, y=mask.grid$y)
		
		values.dif <- result.extract.interpolate(mask.grid=mask.grid, 
					values=values.dif, alt=mask.grid$alt, 
					x=mask.grid$x, y=mask.grid$y)

		if (!silent){ cat("Interpolation done!\n") }
	}

	# create data.frame
	result.values <- data.frame(doy.model=values.model, 
					doy.observed=values.observed, 
					doy.dif=values.dif, 
					x=mask.grid$x, y=mask.grid$y)

	return(result.values)
}
