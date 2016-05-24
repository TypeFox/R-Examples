
if (!isGeneric("erase")) {
	setGeneric("erase", function(x, y, ...)
		standardGeneric("erase"))
}	

.gDif <- function(x, y) {
	xln <- length(x@polygons)
	yln <- length(y@polygons)
	if (xln==0 | yln==0) {
		return(x)
	}
	rn <- row.names(x)
	for (i in xln:1) {
		z <- x[i,]
		for (j in 1:yln) {
			z <- rgeos::gDifference(z, y[j,])
			if (is.null(z)) {
				break
			}
		}
		if (is.null(z)) {
			x <- x[-i,]
			rn <- rn[-i]
		} else {
			x@polygons[i] <- z@polygons
		}
	}
	if (length(x) > 0) {
		j <- which(rgeos::gIsValid(x, byid=TRUE, reason=FALSE))
	#j <- which(gArea(x, byid=TRUE) > 0)			
		x <- x[j,]
		rn <- rn[j]			
	
		if (length(rn) > 0) {
			row.names(x) <- rn
		}
	}
	x
}


setMethod(erase, signature(x='SpatialPolygons', y='SpatialPolygons'),
    function(x, y, ...){ 
	
		requireNamespace("rgeos")

		if (! identical(x@proj4string, y@proj4string) ) {
			warning('non identical CRS')
			y@proj4string <- x@proj4string
		}
		
		if (!.hasSlot(x, 'data')) {
			d <- data.frame(ID=1:length(x@polygons))
			rownames(d) <- row.names(x)
			x <- SpatialPolygonsDataFrame(x, data=d)
			dropframe <- TRUE
		} else {
			dropframe <- FALSE
		}

		y <- aggregate(y)
		
		int <- rgeos::gIntersects(x, y, byid=TRUE)
		int1 <- apply(int, 2, any)
		int2 <- apply(int, 1, any)
				
		if (sum(int1) == 0) { # no intersections
			return(x)
		}
		
		if (all(int1)) {
			part1 <- NULL
		} else {
			part1 <- x[!int1,]
		}
		part2 <- .gDif(x[int1,], y[int2,])

		part2 <- SpatialPolygonsDataFrame(part2, x@data[match(row.names(part2), rownames(x@data)), ,drop=FALSE])
		if (!is.null(part1)) {
			part2 <- rbind(part1, part2)
		}
			
		if (length(part2@polygons) > 1) {	
			part2 <- aggregate(part2, colnames(part2@data))
		}
		if (dropframe) {
			return( as(part2, 'SpatialPolygons') )
		} else {
			return( part2 )
		}
	}
)

