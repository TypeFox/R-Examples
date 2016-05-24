ova2spat <- function( # This function converts an ovariable into a SpatialPointsDataFrame.
		ovariable, # An evaluated ovariable that has coordinate indices.
		coords, # The names of the coordinate indices as a character vector, first x then y.
		proj4string # Projection identifier or specification as character string. See http://spatialreference.org/
) {
	temp <- ovariable@output
	
	# Transform coordinates into numeric format.
	
	for(i in coords) { 
		if(is(temp[[i]], "factor"))    temp[[i]] <- levels(temp[[i]])[temp[[i]]]
		if(is(temp[[i]], "character")) temp[[i]] <- as.numeric(temp[[i]])
	}
	
	# Define the coordinate points first, then add other ovariable output to it.
	
	sp <- SpatialPoints(temp[coords], CRS(proj4string))
	out <- SpatialPointsDataFrame(sp, temp[!colnames(temp) %in% coords])
	
	#Transform the projection to longitude-latitude system.
	
	epsg4326String <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
	out <- spTransform(out,epsg4326String)
	
	return(out)
}