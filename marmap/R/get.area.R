get.area <- function(mat, level.inf, level.sup=0, xlim=NULL, ylim=NULL) {
	
	# require(geosphere)
	
	if (!is(mat,"bathy") & !is(mat,"buffer")) stop("mat must be of class 'bathy' or 'buffer'")
	
	if (is(mat,"buffer")) mat <- mat[[1]]
	
	lon <- as.numeric(rownames(mat))
	lat <- as.numeric(colnames(mat))

	if (is.null(xlim)) xlim <- range(lon)
	if (is.null(ylim)) ylim <- range(lat)
	if (any(is.na(mat))) mat[is.na(mat)] <- 100000*max(mat,na.rm=TRUE)

	x1b <- which.min(abs(lon - xlim[1]))
	y1b <- which.min(abs(lat - ylim[1]))
	x2b <- which.min(abs(lon - xlim[2]))
	y2b <- which.min(abs(lat - ylim[2]))

	mat <- mat[x1b:x2b, y1b:y2b]
	lon <- as.numeric(rownames(mat))
	lat <- as.numeric(colnames(mat))
	
	cell.width <- (lon[2] - lon[1])
	cell.height <- (lat[2] - lat[1])
	
	bathy2 <- ifelse(mat >= level.inf & mat <= level.sup, 1, 0)
	cells <- apply(bathy2,2,sum)
	cells <- cells[cells!=0]
	c.lat <- as.numeric(names(cells))
	
	poly <- list()
	for (i in 1:length(cells))
		poly[[i]] <- rbind(	c(0,c.lat[i]),
							c(0,c.lat[i]+cell.height),
							c(cells[i]*cell.width,c.lat[i]+cell.height),
							c(cells[i]*cell.width,c.lat[i])
							)
	
	surf <- sum(sapply(poly,geosphere::areaPolygon,r=6371))
	return <- list(Square.Km=surf, Area=bathy2, Lon=lon, Lat=lat)
	
}