#'Generate a hemiphoto object
#'
#'@description Construct an object that contains information on shading by the canopy, as
#'measured with hemi-spherical photographs ('hemiphotos'). Reads one of two
#'formats. See YplantQMC website for example files, and the instruction manual
#'for more background.
#'
#'
#'@aliases setHemi evalHemi plot.yphemi
#'@param canfile A canopy file (see Vignette for format example).
#'@param canopytrans Minimum transmission of canopy (when gap fraction = 0).
#'@param hemi An object generated with \code{setHemi}.
#'@param altitude,azimuth Viewing altitude and azimuth.
#'@param met Optionally, a \code{ypmet} object.
#'@param degrees Whether altitude and azimuth were given in degrees (if FALSE,
#'radians).
#'@param x For \code{plot.yphemi}, a \code{yphemi} object
#'@param sungap If TRUE, and a \code{ypmet} object is used, plots the sun along
#'the path, with the size relative to the gap fraction
#'@param projection If 'iso', each altitude bin has the same width in the plot
#'(Default). If 'flat', projection is as seen from above.
#'@param warn If TRUE and sungap=TRUE, warns when gap fractio is very low.
#'@param bordercol Color of the grid separating the sky sectors. Use \code{NA}
#'to omit borders.
#'@param \dots Further arguments passed to plot.default
#'@author Remko Duursma
#'@seealso \code{\link{setMet}}
#'@references For background on hemiphotos :
#'\url{http://en.wikipedia.org/wiki/Hemispherical_photography}
#'@keywords misc
#'@rdname setHemi
#'@export
setHemi <- function(canfile, canopytrans=0.02){
	
	
	if(!is.matrix(canfile)){
		
		extension <- tolower(substr(canfile, nchar(canfile)-2, nchar(canfile)))
		if(extension == "can")
			hemimat <- readCAN(canfile)
		if(extension == "csv"){	
			hemimat <- as.matrix(read.csv(canfile))	
			if(ncol(hemimat) == 1)
				hemimat <- as.matrix(read.table(canfile, header=TRUE))
			if(ncol(hemimat) == 1)
				stop("Error reading 'canfile' - save as comma or space delimited text file.")
			colnames(hemimat) <- NULL
			zenith <- tolower(hemimat[1,1]) == "zenith" # TRUE if zenith given.
			if(!zenith && tolower(hemimat[1,1]) != "altitude")
				stop("Row 1, Column 1 must be 'zenith' or 'altitude'.")
			hemimat <- hemimat[-1,-1]
			storage.mode(hemimat) <- "numeric"
			if(zenith)hemimat <- hemimat[nrow(hemimat):1,]
		}
		if(!extension %in% c("csv","can"))
			stop("Only .CSV or .CAN are currently supported by setHemi")
		
	} else {
		# canfile is already a matrix.
		hemimat <- canfile 
	}	

	# As in Yplant, set the minimum transmission to a parameter,
	# the 'Canopy Transmission Coefficient'.
	hemimat <- pmax(hemimat, canopytrans)
	
	# Number of bins in altitude and azimuth directions.
	nalt <- nrow(hemimat)
	naz <- ncol(hemimat)
	
	# Assume in increasing order; starting at zero (altitude and azimuth);
	# equally spaced bins (angle-wise). Values in 1st col and row NOT read,
	# but should be there anyway (can be empty).
	azbins <- seq(0,2*pi,length=naz+1)
	altbins <- seq(0,pi/2,length=nalt+1)

	hemi <- list()
	hemi$m <- hemimat
	hemi$nalt <- nalt
	hemi$naz <- naz
	hemi$azbins <- azbins
	hemi$altbins <- altbins

	# Set weights for diffuse calculations (i.e., the size of the hemi-tile).
	# Tiles are larger for lower solar altitude.
	w <- c()
	for(i in 1:hemi$nalt)w[i] <- sin(hemi$altbins[i+1]) - sin(hemi$altbins[i])
	weight <- matrix(rep(w,naz),ncol=naz)
	weight <- weight / naz
	
	
	if(!is.matrix(canfile))
		hemi$canfile <- canfile
	else
		hemi$canfile <- NA
	
	class(hemi) <- c("yphemi")
	
return(hemi)
}

