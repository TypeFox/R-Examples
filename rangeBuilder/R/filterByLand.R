# Function to filter occurrences by land/ocean

# coords: 2 column matrix, numeric vector length 2, SpatialPoints object
# returnGood: if TRUE, index of points that pass filter is returned, if FALSE, index of points that fail is returned
# proj: proj4string of input coords. Ignored if input coords are spatial object

# Returns indices of points that make the filter

filterByLand <- function(coords, returnGood = TRUE, proj = '+proj=longlat +datum=WGS84') {

	# if vector, convert to matrix
	if (is.numeric(coords)) {
		coords <- matrix(coords, nrow = 1, ncol = 2)
	}

	if (class(coords) == 'data.frame') {
		coords <- as.matrix(coords)
	}

	if (class(coords) %in% c('SpatialPoints', 'SpatialPointsDataFrame')) {
		
		if (is.na(proj4string(coords))) {
			stop('proj4string must be specified for spatial input.')
		}

		# transform if needed
		if (proj4string(worldRaster) != proj) {
			coords <- sp::spTransform(coords, CRS(proj4string(worldRaster)))
		}

		# convert to matrix
		coords <- as.matrix(as.data.frame(coords))
	}

	if (class(coords) == 'matrix' & proj4string(worldRaster) != proj) {

		#transform
		coords <- SpatialPoints(coords, CRS(proj))
		coords <- sp::spTransform(coords, CRS(proj4string(worldRaster)))
		coords <- as.matrix(as.data.frame(coords))
	}

	if (class(coords) != 'matrix' | mode(coords) != 'numeric') {
		stop('coords must be a numeric matrix.')
	}

	#extract worldRaster values
	e <- extract(worldRaster, coords)
	
	if (returnGood) {
		return(which(!is.na(e)))
	} else {
		return(which(is.na(e)))
	}
}