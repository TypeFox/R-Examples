dist2isobath <- function(mat, x, y=NULL, isobath=0, locator=FALSE, ...) {
	
	if (!is(mat,"bathy")) stop("'mat' must be of class 'bathy'")

	if (!is.numeric(isobath)) stop("'isobath' must be numeric")
	if (length(isobath) >1) {
		isobath <- isobath[1]
		warning("'isobath' must be a single depth or altitude value. Only the first value was used.")
	}

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
		
	coord <- data.frame(x = coord$x, y = coord$y)
		
	# Get contour lines for a given isobath
	lon <- unique(as.numeric(rownames(mat)))
	lat <- unique(as.numeric(colnames(mat)))
	iso <- contourLines(lon, lat, mat, levels = isobath)
		
	# Transform the list from contourLines into a SpatialLines
	iso <- lapply(iso, function(k) sp::Line(matrix(c(k$x,k$y),ncol=2)))
	iso <- sp::SpatialLines(list(sp::Lines(iso, ID = "a")))
	
	# Compute the shortest great circle distance between each point and the isobath
	d <- suppressWarnings(geosphere::dist2Line(coord,iso))
	
	d <- data.frame(d[,1],coord,d[,2:3])
	colnames(d) <- c("distance", "start.lon", "start.lat", "end.lon", "end.lat")
	return(d)
}