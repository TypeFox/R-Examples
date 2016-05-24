readShapePoints <- function(fn, proj4string=CRS(as.character(NA)), 
	verbose=FALSE, repair=FALSE) {
	suppressWarnings(Map <- read.shape(filen=fn, verbose=verbose,
	    repair=repair))
	suppressWarnings(.Map2SPDF(Map, proj4string=proj4string))
}

.Map2SPDF <- function(Map, IDs, proj4string=CRS(as.character(NA))) {
	if (missing(IDs))
		IDs <- as.character(sapply(Map$Shapes, function(x) x$shpID))
	coords <- Map2points(Map)
	oldClass(coords) <- NULL
	rownames(coords) <- IDs
	attr(coords, "shpID") <- NULL
	attr(coords, "maplim") <- NULL

	df <- Map$att.data
	rownames(df) <- IDs
	res <- SpatialPointsDataFrame(coords=coords, data=df, 
		proj4string=proj4string, match.ID=TRUE)
	res
}

writePointsShape <- function(x, fn, factor2char = TRUE, max_nchar=254) {
        stopifnot(is(x, "SpatialPointsDataFrame"))
	df <- as(x, "data.frame")
	coords <- coordinates(x)
	suppressWarnings(write.pointShape(coordinates=coords, df=df, file=fn,  
		factor2char = factor2char, max_nchar=max_nchar))
}


