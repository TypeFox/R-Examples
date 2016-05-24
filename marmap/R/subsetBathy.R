subsetBathy <- function(mat, x, y=NULL, locator=TRUE, ...) {
	
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
			cat('Waiting for interactive input: click any number of times on the map\n')
			coord <- locator(type="o",...)
		}

		as.numeric(rownames(mat)) -> lon
		as.numeric(colnames(mat)) -> lat

		outside.lon <- any(findInterval(coord$x,range(lon),rightmost.closed=TRUE) != 1)
		outside.lat <- any(findInterval(coord$y,range(lat),rightmost.closed=TRUE) != 1)
		if (outside.lon | outside.lat) stop("Some data points are oustide the range of mat")

		out <- data.frame(Lon=coord$x, Lat=coord$y)

		if (nrow(out)==1) stop("'subset.bathy' needs at least two points")

		if (nrow(out)==2) {
			try(rect(min(out$Lon),min(out$Lat),max(out$Lon),max(out$Lat)),silent=TRUE)
			x1 <- which(abs(lon-out$Lon[1])==min(abs(lon-out$Lon[1])))
			y1 <- which(abs(lat-out$Lat[1])==min(abs(lat-out$Lat[1])))
			x2 <- which(abs(lon-out$Lon[2])==min(abs(lon-out$Lon[2])))
			y2 <- which(abs(lat-out$Lat[2])==min(abs(lat-out$Lat[2])))
			new.bathy <- mat[x1:x2, y1:y2]
			new.bathy <- check.bathy(new.bathy)
			class(new.bathy) <- "bathy"
			
		}

		if (nrow(out)>2) {
			x1 <- which(abs(lon-min(out$Lon))==min(abs(lon-min(out$Lon))))
			y1 <- which(abs(lat-min(out$Lat))==min(abs(lat-min(out$Lat))))
			x2 <- which(abs(lon-max(out$Lon))==min(abs(lon-max(out$Lon))))
			y2 <- which(abs(lat-max(out$Lat))==min(abs(lat-max(out$Lat))))
			new.bathy <- mat[x1:x2, y1:y2]
			new.bathy <- check.bathy(new.bathy)
			class(new.bathy) <- "bathy"
		
			xyz <- as.matrix(as.xyz(new.bathy))
			out <- as.matrix(out)
			inside <- sp::point.in.polygon(xyz[,1],xyz[,2],out[,1],out[,2])
			xyz[inside==0,3] <- NA
			new.bathy <- as.bathy(xyz)
		}

	
		return(new.bathy)

}
