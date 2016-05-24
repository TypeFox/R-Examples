get.depth <- function(mat, x, y=NULL, locator=TRUE, distance=FALSE, ...){

	if (!is(mat,"bathy")) stop("'mat' must be of class 'bathy'")

	if (locator == FALSE) {

		if (!is.null(y) & !is.vector(y)) stop("'y' must be a numeric vector or NULL")
		if (!is.null(y) & !is.numeric(y)) stop("'y' must be a numeric vector or NULL")

		if (is.list(x)) {
			if (length(x)!=2) stop("if 'x' is a list, it must contain only two vectors of the same length (longitude and latitude)")
			if (!is.vector(x[[1]]) | !is.vector(x[[2]])) stop("if 'x' is a list, it must contain only two vectors of the same length (longitude and latitude)")
			if (length(x[[1]]) != length(x[[2]])) stop("if 'x' is a list, it must contain only two vectors of the same length (longitude and latitude)")
			if (!is.null(y)) warning("'y' has been ignored\n")
			coord <- x ; names(coord) <- c("x","y")
			}
		
		if (is.data.frame(x) | is.matrix(x)) {
			
			x <- as.matrix(x)

			if (ncol(x) > 2) {
				warning("'x' has too many columns. Only the first two will be considered\n")
				x <- x[,1:2]
				coord <- list(x=x[,1],y=x[,2])
				if (!is.null(y)) warning("only the first two columns of 'x' were considered. 'y' has been ignored\n")
				}

			if (ncol(x) == 2) {
				coord <- list(x=x[,1],y=x[,2])
				if (!is.null(y)) warning("since 'x' has 2 columns, 'y' has been ignored\n")
				}

			if (ncol(x) == 1) {
				if (is.null(y)) stop("with 'locator=FALSE', you must supply both 'x' and 'y' or a 2-column matrix-like 'x'")
				coord <- list(x=x,y=y)
				} 

			}
		
		if (!is.list(x)) {	
			if (is.vector(x) & !is.numeric(x)) stop("'x' must be numeric")
			if (is.vector(x) & is.numeric(x)) {
				if (is.null(y)) stop("with 'locator=FALSE', you must either provide both 'x' and 'y' or a 2-column matrix-like 'x'")
					if (length(x) != length(y)) warning("The lengths of 'x' and 'y' differ. Some values have been recycled\n")
						coord <- list(x=x,y=y)
				}
			}
			
		} else {
			cat("Waiting for interactive input: click any number of times on the map, then press 'Esc'\n")
			coord <- locator(type="p",...)
		}
	
	as.numeric(rownames(mat)) -> lon
	as.numeric(colnames(mat)) -> lat
	
	outside.lon <- any(findInterval(coord$x,range(lon),rightmost.closed=TRUE) != 1)
	outside.lat <- any(findInterval(coord$y,range(lat),rightmost.closed=TRUE) != 1)
	
	if (outside.lon | outside.lat) stop("Some data points are oustide the range of mat")
		
	out <- data.frame(lon=coord$x, lat=coord$y)
	out$depth <- apply(out, 1, function(x) mat[ which(abs(lon-x[1])==min(abs(lon-x[1]))) , which(abs(lat-x[2])==min(abs(lat-x[2]))) ][1])
	
	if(distance){
		
		if (nrow(out) == 1) stop("Cannot compute distance with only one point. Either set distance=FALSE or add more points")
		
		deg2km <- function(x1, y1, x2, y2) {

			x1 <- x1*pi/180
			y1 <- y1*pi/180
			x2 <- x2*pi/180
			y2 <- y2*pi/180

			dx <- x2-x1
			dy <- y2-y1

			fo <- sin(dy/2)^2 + cos(y1) * cos(y2) * sin(dx/2)^2
			fos <- 2 * asin(min(1,sqrt(fo)))

			return(6371 * fos)
		}
		
		dist.km = NULL
		for(i in 2:length(out$depth)){
			dist.km = c(dist.km, deg2km(x1=out$lon[i-1],y1=out$lat[i-1],x2=out$lon[i],y2=out$lat[i]))
		}
		out$dist.km <- cumsum(c(0,dist.km))
		out <- out[,c(1,2,4,3)]
	}
	
	return(out)

}
