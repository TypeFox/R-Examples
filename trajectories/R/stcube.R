# Provide stcube methods.

map3d = function(map, z, ...) {
	if (!requireNamespace("rgl", quietly = TRUE))
		stop("rgl required")
	if(length(map$tiles) != 1)
		stop("Pass single map tile only.")
	nx = map$tiles[[1]]$xres
	ny = map$tiles[[1]]$yres
	xmin = map$tiles[[1]]$bbox$p1[1]
	xmax = map$tiles[[1]]$bbox$p2[1]
	ymin = map$tiles[[1]]$bbox$p1[2]
	ymax = map$tiles[[1]]$bbox$p2[2]
	xc = seq(xmin, xmax, len = ny)
	yc = seq(ymin, ymax, len = nx)
	col = matrix(data = map$tiles[[1]]$colorData, nrow = ny, ncol = nx)
	m = matrix(data = z, nrow = ny, ncol = nx)
	
	rgl::surface3d(x = xc, y = yc, z = m, col = col, lit = FALSE, ...)
}

normalize = function(time, by = "week") {
	tn = as.numeric(time)

	switch(by,
		minute = (tn %% 60),
		hour = (tn %% 3600) / 60 , # decimal minute of the hour
		day = (tn %% (3600 * 24)) / 3600, # decimal hour of the day
		week = (tn %% (3600 * 24 * 7)) / 24 / 3600, # decimal day of the week
		stop(paste("unknown value for by: ", by)))
}

OSM = function(xlim, ylim, mapZoom, mapType, projection) {
	if (!requireNamespace("OpenStreetMap", quietly = TRUE))
		stop("package OpenStreetMap required")
	map = OpenStreetMap::openmap(upperLeft = c(ylim[2], xlim[1]),
		lowerRight = c(ylim[1], xlim[2]), zoom = mapZoom, type = mapType)
	OpenStreetMap::openproj(x = map, projection = projection)
}

if(!isGeneric("stcube"))
	setGeneric("stcube", function(x, ...)
		standardGeneric("stcube"))

setMethod("stcube", signature(x = "Track"),
	function(x, xlab = "x", ylab = "y", zlab = "t", type = "l", aspect, 
		xlim = stbox(x)[[1]], ylim = stbox(x)[[2]], zlim = stbox(x)$time, 
		showMap = FALSE, mapType = "osm", mapZoom = NULL, ..., y, z) {
		# "y" and "z" are ignored, but added to avoid ... absorbs them
		if (!requireNamespace("rgl", quietly = TRUE))
			stop("rgl required")
		coords = coordinates(x@sp)
		time = index(x@time)
		time <- time - min(time) # seconds from start
		if(missing(aspect))
			aspect = if((asp = mapasp(x@sp)) == "iso") "iso" else c(1, asp, 1)
		if (missing(zlim))
			zlim = range(time)
		rgl::plot3d(x = coords[, 1], y = coords[, 2], z = time, xlab = xlab,
			ylab = ylab, zlab = zlab, type = type, aspect = aspect, xlim = xlim,
			ylim = ylim, zlim = zlim, ...)
		if(showMap)
			map3d(map = OSM(xlim, ylim, mapZoom, mapType, proj4string(x)), z = time[1])
	}
)

setMethod("stcube", signature(x = "Tracks"),
	function(x, xlab = "x", ylab = "y", zlab = "t", type = "l", aspect,
		xlim = stbox(x)[[1]], ylim = stbox(x)[[2]], zlim = stbox(x)$time, 
		showMap = FALSE, mapType = "osm", normalizeBy = "week",
		mapZoom = NULL, ..., y, z, col) {
		# "y" and "z" are ignored, but added to avoid ... absorbs them
		if (!requireNamespace("rgl", quietly = TRUE))
			stop("rgl required")
		dim = dim(x@tracks[[1]])["geometries"]
		coordsAll = do.call(rbind, lapply(x@tracks, function(x) coordinates(x@sp)))
		timeAll = normalize(do.call(c, lapply(x@tracks,
			function(x) index(x@time))), normalizeBy)
		col = rainbow(length(x@tracks))
		if(missing(aspect))
			# mapasp() processes objects of class Spatial* only.
			aspect = if((asp = mapasp(as(x, "SpatialLines"))) == "iso") "iso" else c(1, asp, 1)
		if (missing(zlim))
			zlim = range(timeAll)
		rgl::plot3d(x = coordsAll[1:dim, 1], y = coordsAll[1:dim, 2],
			z = timeAll[1:dim], xlab = xlab, ylab = ylab, zlab = zlab,
			type = type, col = col[1], aspect = aspect, xlim = xlim,
			ylim = ylim, zlim = zlim, ...)
		tracks = x@tracks[-1]
		for(t in seq_along(tracks)) {
			coords = coordinates(tracks[[t]]@sp)
			time = normalize(index(tracks[[t]]@time), normalizeBy)
			rgl::lines3d(x = coords[, 1], y = coords[, 2], z = time, col = col[t+1])
		}
		if(showMap)
			map3d(map = OSM(xlim, ylim, mapZoom, mapType, proj4string(x)), z = timeAll[1])
	}
)

setMethod("stcube", signature(x = "TracksCollection"),
	function(x, xlab = "x", ylab = "y", zlab = "t", type = "l", aspect,
		xlim = stbox(x)[[1]], ylim = stbox(x)[[2]], zlim = stbox(x)$time, 
		showMap = FALSE, mapType = "osm", normalizeBy = "week",
		mapZoom = NULL, ..., y, z, col) {
		# "y", "z" and "col" are ignored, but included in the method signature
		if (!requireNamespace("rgl", quietly = TRUE))
			stop("rgl required")
		dim = dim(x@tracksCollection[[1]]@tracks[[1]])["geometries"]
		coordsAll = do.call(rbind, lapply(x@tracksCollection,
			function(x) do.call(rbind, lapply(x@tracks, function(y) coordinates(y@sp)))))
		timeAll = normalize(do.call(c, lapply(x@tracksCollection,
			function(x) do.call(c, lapply(x@tracks,
				function(y) index(y@time))))), normalizeBy)
		if (missing(col))
			col = rainbow(length(x@tracksCollection))
		if(missing(aspect))
			# mapasp() processes objects of class Spatial* only.
			aspect = if((asp = mapasp(as(x, "SpatialLines"))) == "iso") "iso" else c(1, asp, 1)
		if (missing(zlim))
			zlim = range(timeAll)
		rgl::plot3d(x = coordsAll[1:dim, 1], y = coordsAll[1:dim, 2],
			z = timeAll[1:dim], xlab = xlab, ylab = ylab, zlab = zlab,
			type = type, col = col[1], aspect = aspect, xlim = xlim,
			ylim = ylim, zlim = zlim, ...)
		for(tz in seq_along(x@tracksCollection)) {
			if(tz == 1)
				tracks = x@tracksCollection[[tz]]@tracks[-1]
			else
				tracks = x@tracksCollection[[tz]]@tracks
			for(t in seq_along(tracks)) {
				coords = coordinates(tracks[[t]]@sp)
				time = normalize(index(tracks[[t]]@time), normalizeBy)
				rgl::lines3d(x = coords[, 1], y = coords[, 2], z = time, col = col[tz])
			}
		}
		if(showMap)
			map3d(map = OSM(xlim, ylim, mapZoom, mapType, proj4string(x)), z = timeAll[1])
	}
)
