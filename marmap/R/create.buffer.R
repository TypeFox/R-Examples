create.buffer <- function(x, loc, radius, km=FALSE){
	
	mat <- x # S3 compatibility
	
	if (!is(mat,"bathy")) stop("x must be an object of class bathy")
	if (!is.data.frame(loc)) stop("loc must be a two-column data.frame (longitude and latitude)")
	if (!is.numeric(radius)) stop("radius must be numeric")
	if (length(radius) > 1) warning("only the first value of radius was used")
	
	xyz <- as.xyz(mat)
	
	if (km) {
		radius.km <- radius
		radius <- 180 * radius.km/(pi*6372.798)
	} else {
		radius.km <- radius*pi*6372.798/180
	}
	
	map <- sp::SpatialPixelsDataFrame(points = xyz[,1:2], data = xyz, tolerance = 0.006)
    lo2 <- sp::SpatialPointsDataFrame(loc, data = loc)
	
	adehabitatMA::adeoptions(epsilon=0.001)
	temp <- adehabitatMA::buffer(lo2, map, radius)
	adehabitatMA::adeoptions(epsilon=0.00000001)
		
	temp <- -as.bathy(as(temp,'SpatialGridDataFrame'))
	
	mat[temp==0] <- NA

	out <- list(buffer = mat, center = loc, radius = radius, radius.km = radius.km)
	class(out) <- "buffer"
	return(out)
	
}