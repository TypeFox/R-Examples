# function to determine overall extent of a list of polygons
getExtentOfList <- function(shapes) {	
	x <- lapply(shapes, function(x) bbox(x))
	minLong <- min(sapply(x, function(x) x[1], simplify = TRUE))
	maxLong <- max(sapply(x, function(x) x[3], simplify = TRUE))
	minLat <- min(sapply(x, function(x) x[2], simplify = TRUE))
	maxLat <- max(sapply(x, function(x) x[4], simplify = TRUE))
	return(list(minLong = minLong, maxLong = maxLong, minLat = minLat, maxLat = maxLat))
}	
