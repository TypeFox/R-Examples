#'@rdname setHemi
evalHemi <- function(hemi, altitude=NULL, azimuth=NULL, 
	met=NULL, degrees=TRUE){

	if(!inherits(hemi, "yphemi"))
		stop("'hemi' needs to be an YplantQMC hemiphoto object (see ?setHemi).")

	if(is.null(met)){
		if(degrees){
			altout <- altitude * pi/180
			azout <- azimuth * pi/180
		} else {
			altout <- altitude
			azout <- azimuth
		}
		if(length(altout) != length(azout))
			stop("Altitude (altout) and azimuth (azout) need to be vectors of equal length.")
	} else {  # met not NULL
		altout <- pi/180 * met$dat$altitude
		azout <- pi/180 * met$dat$azimuth
	}

	# Look up gapfraction for the hemiphoto tile.
	p_azbin <- findInterval(azout, hemi$azbins)
	p_altbin <- findInterval(altout, hemi$altbins)
  
	# Sometimes a point is numerically > pi/2; set altitude to max bin.
	p_altbin[p_altbin > hemi$nalt] <- hemi$nalt
	
  # Sometimes azimuth is exactly 360 (2*pi), ends in last bin. Put in first.
	p_azbin[p_azbin > hemi$naz] <- 1
	
	intgapfrac <- c()
	for(i in 1:length(p_azbin))
		intgapfrac[i] <- hemi$m[p_altbin[i], p_azbin[i]]
		
return(data.frame(altitude=altout, azimuth=azout, gapfraction=intgapfrac))
}	
	

	
	
	
	
	