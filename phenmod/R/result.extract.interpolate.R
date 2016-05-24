result.extract.interpolate <- function(mask.grid, values, alt, x, y){
	# remove invalid values
	valid.values <- which((is.na(values)==FALSE) & (values > -9999))
	values <- values[valid.values]
	alt <- alt[valid.values]
	x <- x[valid.values]
	y <- y[valid.values]

	if (length(valid.values) <= 1){
		warning("Interpolation: not enough valid values")
		return(rep(NA, length(mask.grid)))
	}

	# create data.frame with valid values
	final.data <- data.frame(x=x, y=y, alt=alt, values=values)
	
	# remove old variables
	rm(x,y,alt,values,valid.values)

	# find observations with identical coordinates and modify them
	dubli.x <- duplicated(final.data$x)
	dubli.y <- duplicated(final.data$y)
  	uni <- unique(final.data$x)
	dif <- length(dubli.x) - length(uni)
  	final.data$x[dubli.x] <- round(rnorm(dif, final.data$x[dubli.x], 2), digits=0)
	rm(dubli.x, dubli.y, uni, dif)

	# variogram estimations
	vgm1 <- variogram(values~alt, locations=~x+y, 
			data=final.data, nmax=20, nmin=10)
 	model.1 <- fit.variogram(vgm1,vgm(1000,"Sph",200000,600))

	# external drift kriging
	values.edk <-  krige(values~alt, locations=~x+y, 
				data=final.data, newdata=mask.grid, model=model.1, 
				nmax=20, nmin=10, maxdist=200000, na.action=na.pass)

	values.edk <- values.edk$var1.pred

	return(values.edk)
}