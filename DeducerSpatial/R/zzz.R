# TODO: Add comment
# 
# Author: ianfellows
###############################################################################

.containedBy <- function (minLat, minLon, maxLat, maxLon, coords) {
	minMerc <- projectMercator(minLat, minLon)
	maxMerc <- projectMercator(maxLat, maxLon)
	
	#print(minMerc)
	#print(maxMerc)
	
	minLat <- minMerc[[1]]
	maxLat <- maxMerc[[1]]
	
	minLon <- minMerc[[2]]
	maxLon <- maxMerc[[2]]
	
	for (i in 1:dim(coords)[[1]]) { 
		lat <- coords[i, 1] 
		lon <- coords[i, 2]
		if (!( lat >= minLat && lat <= maxLat && lon >= minLon && lon <= maxLon)) {
			return(FALSE)
		}
	}
	return(TRUE)
}

#Used for rectangle subsetting
#Returns true if *ALL* of the coords are contained by the rectangle
.containedBy2 <- function (minMerc, maxMerc, coords) {
	minLat <- minMerc[[1]]
	maxLat <- maxMerc[[1]]
	minLon <- minMerc[[2]]
	maxLon <- maxMerc[[2]]
	all(coords[,1] >= minLat & coords[,1] <= maxLat & coords[,2] >= minLon & coords[,2] <= maxLon)
}

#For Polygons
#TODO handle error when no polys left.  Probably should do nothing.
#TODO These following functions could likely be consolidated
.subsetPoly <- function (minLat, minLon, maxLat, maxLon, polyDf, removeSelection) {
	minMerc <- projectMercator(minLat, minLon)
	maxMerc <- projectMercator(maxLat, maxLon)
	# The XOR inverts the function in any easy way w/o if/else statements
	.contained <- function(poly){return(xor(removeSelection, .containedBy2(minMerc, maxMerc, poly)))}
	
	nr <- nrow(polyDf)
	indices <- rep(FALSE,nr)
	
	for (i in 1:nr){
		indices[i] <- .contained(t(polyDf[i,]@bbox)) 
	}
	
	if(all(!indices))
		return(NULL)
	polyDf[indices,]
}

.subsetLines <- function (minLat, minLon, maxLat, maxLon, polyDf, removeSelection) {
	minMerc <- projectMercator(minLat, minLon)
	maxMerc <- projectMercator(maxLat, maxLon)
	# The XOR inverts the function in any easy way w/o if/else statements
	.contained <- function(poly){return(xor(removeSelection, .containedBy2(minMerc, maxMerc, poly@coords)))}
	
	dupDf <- polyDf
	
	for (i in 1:length(polyDf@lines)){# each list of polygons
		dupDf@lines[[i]]@Lines <- Filter(.contained, polyDf@lines[[i]]@Lines)
		#.containedBy(minLat, minLon, maxLat, maxLon, poly@coords)
	}
	
	indices <- 1:length(polyDf)
	indices <- Filter(function(x){return(length(dupDf@lines[[x]]@Lines) > 0)}, indices)
	if (length(indices) == 0)
	{
		return(NULL)
	}
	else
	{
		dupDf <- polyDf[indices, ]
		return(dupDf)
	}
}

.subsetPoints <- function (minLat, minLon, maxLat, maxLon, pointsDf, removeSelection) {
	minMerc <- projectMercator(minLat, minLon)
	maxMerc <- projectMercator(maxLat, maxLon)
	dupDf <- pointsDf
	
	.contained <- function(x) {
		return(xor(removeSelection, .containedBy2(minMerc, maxMerc, dupDf[x,]@coords) ) )
	}
	
	indices <- 1:length(pointsDf)
	indices <- Filter(.contained, indices)
	
	if (length(indices) == 0)
	{
		return(NULL)
	}
	else
	{
		dupDf <- pointsDf[indices, ]
		return(dupDf)	
	}
}



.onLoad <- function(libname, pkgname) {

	deducerLoaded <- try(.deducer != .jnull(),silent=TRUE)
	if(inherits(deducerLoaded,"try-error") || !deducerLoaded)
		return(NULL)
	
	if (!nzchar(Sys.getenv("NOAWT")) || .jgr==TRUE){
		.jpackage(pkgname,lib.loc=libname)  
		.jengine(TRUE)
		DeducerSpatial <- J("edu.cens.spatial.DeducerSpatial")
		DeducerSpatial$init()
		.registerDialog("Load Census Data",.makeCensusDialog)
	}
	
	#x <- .jnew(J("edu.cens.spatial.Spatial"));
	#x <- try(x$initJGR(), silent=TRUE);
	return(TRUE)
}

