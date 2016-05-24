# Segments are data.frames with a segment on each row, with x0 y0 x1 y1 the
# first four values, followed by attributes.

setClass("segments", contains = "data.frame")

# Coerce to segments.

setAs("Track", "segments", 
	function(from) {
		cc = coordinates(from@sp)
		t = index(from@time)
		df = from@data
		data.frame(x0 = head(cc[,1], -1), y0 = head(cc[,2], -1),
			x1 = tail(cc[,1], -1), y1 = tail(cc[,2], -1),
			time = head(t, -1), head(df, -1), from@connections)
	}
)

setAs("Tracks", "segments", function(from) {
		ret = do.call(rbind, lapply(from@tracks, 
			function(x) as(x, "segments")))
		ret$Track = rep(names(from@tracks), 
			times = sapply(from@tracks, length) - 1)
		ret
	}
)

setAs("TracksCollection", "segments",
	function(from) {
		l = lapply(from@tracksCollection, function(x) as(x, "segments"))
		ret = do.call(rbind, l)
		ret$IDs = rep(names(from@tracksCollection), times = sapply(l, nrow))
		ret
	}
)

# Coerce to data frame.

setAs("Track", "data.frame", 
	function(from) as(as(from, "STIDF"), "data.frame")
)

setAs("Tracks", "data.frame", 
	function(from) {
		l = lapply(from@tracks, function(x) rbind(as(x, "data.frame"), NA))
		d = do.call(rbind, l)
		d$Track = rep(names(from@tracks), times = sapply(l, nrow))
		d
	}
)

setAs("TracksCollection", "data.frame", 
	function(from) {
		l = lapply(from@tracksCollection, function(x) as(x, "data.frame"))
		ret = do.call(rbind, l)
		data.frame(ret, IDs = rep(names(from@tracksCollection), times = sapply(l, nrow)))
	}
)

# Coerce to Line, Lines, SpatialLines and SpatialLinesDataFrame. 

setAs("Track", "Line", 
	function(from) Line(coordinates(from@sp))
)

setAs("Track", "Lines",
	function(from) Lines(list(as(from, "Line")), "ID")
)

setAs("Track", "SpatialLines",
	function(from) SpatialLines(list(as(from, "Lines")), CRS(proj4string(from)))
)

setAs("Tracks", "Lines", 
	function(from) {
		tz = from@tracks
		# The Lines ID is made up of the conjunction of the first and last Track
		# ID, using hyphen as separator.
		Lines(lapply(tz, function(x) as(x, "Line")), 
			paste(names(tz)[1], names(tz)[length(tz)], sep = "-"))
	}
)

setAsWithID <- function(from, ID=NA) {
	l = lapply(from@tracks, function(x) as(x, "Lines"))
	for (i in seq_along(l))
		l[[i]]@ID = paste(switch(is.na(ID),"ID",ID), "_", i, sep="")
	SpatialLines(l, CRS(proj4string(from)))
}

setAs("Tracks", "SpatialLines", function(from) setAsWithID(from))

setAs("Tracks", "SpatialLinesDataFrame", 
	function(from) SpatialLinesDataFrame(as(from, "SpatialLines"), 
		from@tracksData, match.ID = FALSE)
)

setAs("TracksCollection", "SpatialLines", 
	function(from) {
		l <- lapply(from@tracksCollection, function(tracksObj) as(tracksObj, "Lines"))

		if (is.null(rownames(from@tracksCollectionData))) 
			tracksIDs <- paste("ID", 1:length(from@tracksCollection), sep="")
		else 
		tracksIDs <- rownames(from@tracksCollectionData)
		
		trackIDs <- rep(tracksIDs, sapply(from@tracksCollection, length))
		
		for (i in seq_along(l)) {
			l[[i]]@ID <- paste(trackIDs[i], l[[i]]@ID, sep="_")
		}
		
		SpatialLines(l, CRS(proj4string(from)))
	}
)


setAs("TracksCollection", "SpatialLinesDataFrame", 
	function(from) SpatialLinesDataFrame(as(from, "SpatialLines"),
		from@tracksCollectionData, match.ID = FALSE)
)

# Coerce to xts; Track is automatic through STIDF.

setAs("Tracks", "xts",
	function(from)
		do.call(rbind, lapply(from@tracks, function(x) as(x, "xts")))
)

setAs("TracksCollection", "xts",
	function(from)
		do.call(rbind, lapply(from@tracksCollection, function(x) as(x, "xts")))
)

# Coerce to STIDF.

setAs("Tracks", "STIDF",
	function(from)
		do.call(rbind, lapply(from@tracks, function(x) as(x, "STIDF")))
)

setAs("TracksCollection", "STIDF",
	function(from)
		do.call(rbind, lapply(from@tracksCollection, 
			function(x) as(x, "STIDF")))
)

# Coerce to Spatial*

setAs("Track", "Spatial",
	function(from) {
		from@data$time = index(from@time)
		if (!all(from@data$time == from@endTime))
			from@data$endTime = from@endTime
		addAttrToGeom(from@sp, from@data, match.ID = FALSE)
	}
)

setAs("Tracks", "Spatial",
	function(from) {
		ret = do.call(rbind, lapply(from@tracks, 
			function(x) as(x, "Spatial")))
		ret$Track = rep(names(from@tracks), times = lapply(from@tracks, length))
		ret
	}
)

setAs("TracksCollection", "Spatial",
	function(from)
		do.call(rbind, lapply(from@tracksCollection, 
			function(x) as(x, "Spatial")))
)

# Coerce to SpatialPointsDataFrame.

setAs("Track", "SpatialPointsDataFrame", function(from) as(from, "Spatial"))
setAs("Tracks", "SpatialPointsDataFrame", function(from) as(from, "Spatial"))
setAs("TracksCollection", "SpatialPointsDataFrame", function(from) as(from, "Spatial"))

# Provide coordinates methods.

setMethod("coordinates", "Track",
	function(obj) coordinates(obj@sp)
)

setMethod("coordinates", "Tracks",
	function(obj) do.call(rbind, lapply(obj@tracks,
		function(x) coordinates(x)))
)

setMethod("coordinates", "TracksCollection",
	function(obj) do.call(rbind, lapply(obj@tracksCollection,
		function(x) coordinates(x)))
)

# Provide proj4string methods.

setMethod("proj4string", signature(obj = "Track"),
	function(obj) proj4string(obj@sp)
)
setMethod("proj4string", signature(obj = "Tracks"),
	function(obj) proj4string(obj@tracks[[1]])
)
setMethod("proj4string", signature(obj = "TracksCollection"),
	function(obj) proj4string(obj@tracksCollection[[1]])
)

# Provide plot methods. TODO Make more generic.

plot.TracksCollection <- function(x, y, ..., type = 'l', xlim = stbox(x)[,1],
		ylim = stbox(x)[,2], col = 1, lwd = 1, lty =
		1, axes = TRUE, Arrows = FALSE, Segments = FALSE, add = FALSE) {
	sp = x@tracksCollection[[1]]@tracks[[1]]@sp
	# Submitting arrows() and lines() formals to plot prompts warnings, 
	# so they are absorbed here by localPlot() formals.
	localPlot <- function (..., length, angle, code, xpd, lend, ljoin, lmitre)
		plot (...)
	if (! add)
		localPlot(as(sp, "Spatial"), xlim = xlim, ylim = ylim, axes = axes, ...)
	if (axes == FALSE)
		box()
	if (Arrows || Segments) {
		df = as(x, "segments")
		args = list(x0 = df$x0, y0 = df$y0, x1 = df$x1, y1 = df$y1, 
			col = col, lwd = lwd, lty = lty, ...)
		if (Arrows)
			do.call(arrows, args)
		else
			do.call(segments, args)
	} else {
		linesTrack = function(x, y, ...) lines(as(x@sp, "SpatialLines"), ...)
		linesTracks = function(x, y, ...) mapply(linesTrack, x@tracks, ...)
		invisible(mapply(linesTracks, 
			x@tracksCollection, col = col, lwd = lwd, lty = lty, type = type, ...))
	}
}
setMethod("plot", "TracksCollection", plot.TracksCollection)

setMethod("plot", "Tracks", function(x, ...) plot(TracksCollection(list(x)), ...))

setMethod("plot", c("Track", "missing"), function(x, ...) plot(Tracks(list(x)), ...))

# Provide coordnames methods.

setMethod("coordnames", "Track", function(x) coordnames(x@sp))

setMethod("coordnames", "Tracks", function(x) coordnames(x@tracks[[1]]))

setMethod("coordnames", "TracksCollection", function(x) coordnames(x@tracksCollection[[1]]))

# Provide stbox methods.

# AUTOMATIC: fallback to ST (see ST-methods.R)
#setMethod("stbox", "Track",
#	function(obj) {
#		bb = data.frame(t(bbox(obj@sp)))
#		bb$time = range(index(obj@time))
#		rownames(bb) = c("min", "max")
#		bb
#	}
#)

setMethod("stbox", "Tracks",
	function(obj) {
		df = obj@tracksData
		xr = c(min(df$xmin), max(df$xmax))
		yr = c(min(df$ymin), max(df$ymax))
		tr = c(min(df$tmin), max(df$tmax))
		cn = coordnames(obj)
		ret = data.frame(xr, yr, time = tr)
		colnames(ret)[1:2] = cn
		rownames(ret) = c("min", "max")
		ret
	}
)

setMethod("stbox", "TracksCollection",
	function(obj) {
		df = obj@tracksCollectionData
		xr = c(min(df$xmin), max(df$xmax))
		yr = c(min(df$ymin), max(df$ymax))
		tr = c(min(df$tmin), max(df$tmax))
		cn = coordnames(obj)
		ret = data.frame(xr, yr, time = tr)
		colnames(ret)[1:2] = cn
		rownames(ret) = c("min", "max")
		ret
	}
)

# Provide bbox methods.

setMethod("bbox", "Tracks", function(obj) t(stbox(obj)[1:2]))

setMethod("bbox", "TracksCollection", function(obj) t(stbox(obj)[1:2]))

# Provide over methods.

setMethod("over", c("Track", "Spatial"), function(x, y, ...) over(as(x, "SpatialLines"), y, ...))
setMethod("over", c("Tracks", "Spatial"), function(x, y, ...) over(as(x, "SpatialLines"), y, ...))
setMethod("over", c("TracksCollection", "Spatial"), function(x, y, ...) over(as(x, "SpatialLines"), y, ...))

setMethod("over", c("Track", "xts"), function(x, y, ...) over(as(x, "xts"), y, ...))
setMethod("over", c("Tracks", "xts"), function(x, y, ...) over(as(x, "xts"), y, ...))
setMethod("over", c("TracksCollection", "xts"), function(x, y, ...) over(as(x, "xts"), y, ...))

# Provide aggregate methods.

setMethod("aggregate", "Track", 
	function(x, ...) {
		aggregate(as(x, "SpatialPointsDataFrame"), ...)
	}
)

setMethod("aggregate", "Tracks", 
	function(x, ...) {
		aggregate(as(x, "SpatialPointsDataFrame"), ...)
	}
)

setMethod("aggregate", "TracksCollection",
	function(x, ...) {
		aggregate(as(x, "SpatialPointsDataFrame"), ...)
	}
)

# Provide dimension methods.

dim.Track = function(x) c(geometries = length(x@sp))

dim.Tracks = function(x) c(tracks = length(x@tracks),
	geometries = sum(sapply(x@tracks, dim)))

dim.TracksCollection = function(x) c(IDs = length(x@tracksCollection),
	apply(sapply(x@tracksCollection,dim), 1, sum))

# Provide summary methods.

summary.Track = function(object, ...) {
	obj = list()
	obj$class = class(object)
	obj$dim = dim(object)
	obj$stbox = stbox(object)
	obj$sp = summary(object@sp)
	obj$time = summary(object@time)
	obj$data = summary(object@data)
	obj$connections = summary(object@connections)
	class(obj) = "summary.Track"
	obj
}

setMethod("summary", "Track", summary.Track)

print.summary.Track = function(x, ...) {
	cat(paste("Object of class ", x$class, "\n", sep = ""))
	cat(" with Dimensions: (")
	cat(paste(x$dim, collapse = ", "))
	cat(")\n")
	cat("[[stbox]]\n")
	print(x$stbox)
	cat("[[Spatial:]]\n")
	print(x$sp)
	cat("[[Temporal:]]\n")
	print(x$time)
	cat("[[Data attributes:]]\n")
	print(x$data)
	cat("[[Connections:]]\n")
	print(x$connections)
	invisible(x)
}

summary.Tracks = function(object, ...) {
	obj = list()
	obj$class = class(object)
	obj$dim = dim(object)
	obj$stbox = stbox(object)
	obj$sp = summary(do.call(rbind, lapply(object@tracks, function(x) x@sp)))
	obj$time = summary(do.call(rbind, lapply(object@tracks, function(x) x@time)))
	obj$data = summary(do.call(rbind, lapply(object@tracks, function(x) x@data)))
	obj$connections = summary(do.call(rbind, lapply(object@tracks, function(x) x@connections)))
	class(obj) = "summary.Tracks"
	obj
}

setMethod("summary", "Tracks", summary.Tracks)

print.summary.Tracks = function(x, ...) {
	cat(paste("Object of class ", x$class, "\n", sep = ""))
	cat(" with Dimensions (tracks, geometries): (")
	cat(paste(x$dim, collapse = ", "))
	cat(")\n")
	cat("[[stbox]]\n")
	print(x$stbox)
	cat("[[Spatial:]]\n")
	print(x$sp)
	cat("[[Temporal:]]\n")
	print(x$time)
	cat("[[Data attributes:]]\n")
	print(x$data)
	cat("[[Connections:]]\n")
	print(x$connections)
	invisible(x)
}

summary.TracksCollection = function(object, ...) {
	obj = list()
	obj$class = class(object)
	obj$dim = dim(object)
	obj$stbox = stbox(object)
	obj$sp = summary(do.call(rbind, lapply(object@tracksCollection,
		function(x) do.call(rbind, lapply(x@tracks, function(y) y@sp)))))
	obj$time = summary(do.call(rbind, lapply(object@tracksCollection,
		function(x) do.call(rbind, lapply(x@tracks, function(y) y@time)))))
	obj$data = summary(do.call(rbind, lapply(object@tracksCollection,
		function(x) do.call(rbind, lapply(x@tracks, function(y) y@data)))))
	obj$connections = summary(do.call(rbind, lapply(object@tracksCollection,
		function(x) do.call(rbind, lapply(x@tracks, function(y) y@connections)))))
	class(obj) = "summary.TracksCollection"
	obj
}

setMethod("summary", "TracksCollection", summary.TracksCollection)

print.summary.TracksCollection = function(x, ...) {
	cat(paste("Object of class ", x$class, "\n", sep = ""))
	cat(" with Dimensions (IDs, tracks, geometries): (")
	cat(paste(x$dim, collapse = ", "))
	cat(")\n")
	cat("[[stbox]]\n")
	print(x$stbox)
	cat("[[Spatial:]]\n")
	print(x$sp)
	cat("[[Temporal:]]\n")
	print(x$time)
	cat("[[Data attributes:]]\n")
	print(x$data)
	cat("[[Connections:]]\n")
	print(x$connections)
	invisible(x)
}

# Provide selection methods.
subs.Track <- function(x, i, j, ..., drop = TRUE) {
	Track(as(x, "STIDF")[i, j, ..., drop = drop])
}

setMethod("[", "Track", subs.Track)

subs.Tracks <- function(x, i, j, ... , drop = TRUE) {
	if (missing(i))
		i = 1:length(x@tracks)
	else if (is(i, "Spatial"))
		i = which(!is.na(over(x, geometry(i))))
	else if (is.character(i) && length(i) == 1 && !(i %in% names(x@tracks))) # i time:
		i = sapply(x@tracks, function(x) nrow(as(x, "xts")[i]) > 0)
	else if (is.logical(i))
		i = which(i)
	if (length(i) == 1 && drop)
		x@tracks[[i]]
	else {
		if (!any(i))
			NULL
		else 
			Tracks(x@tracks[i], x@tracksData[i, j, drop=FALSE])
	}
}

setMethod("[", "Tracks", subs.Tracks)

subs.TracksCollection <- function(x, i, j, ... , drop = TRUE) {
	if (!missing(j) && is.character(j)) {
		for(tz in seq_along(x@tracksCollection)) {
			for(t in seq_along(x[tz]@tracks)) {
				data = x[tz][t]@data
				connections = x[tz][t]@connections
				if(j %in% names(data))
					data = data[j]
				else
					# An empty data slot is returned if the passed attribute
					# does not exist. The same applies to the connections slot.
					data = data.frame(matrix(nrow = dim(x[tz][t])["geometries"], ncol = 0))
				if(j %in% names(connections))
					connections = connections[j]
				else
					connections = data.frame(matrix(nrow = dim(x[tz][t])["geometries"] - 1, ncol = 0))
				# Write back the just processed data and connection slots.
				x@tracksCollection[[tz]]@tracks[[t]]@data = data
				x@tracksCollection[[tz]]@tracks[[t]]@connections = connections
			}
		}
	}
	if (missing(i))
		s = 1:length(x@tracksCollection)
	else if (is(i, "Spatial"))
		s = which(!is.na(over(x, geometry(i))))
	else if (is.character(i) && length(i) == 1 && !(i %in% names(x@tracksCollection))) # i time:
		i = lapply(x@tracksCollection, function(x) 
			which(sapply(x@tracks, function(x) nrow(as(x, "xts")[i]) > 0)))
	else if (is.logical(i))
		s = which(i)
	else
		s = i
	if (!missing(i) && is.list(i)) { # might have been created by the lappty() above
		stopifnot(all(sapply(i, function(x) is.numeric(x))))
		s = which(sapply(i, function(x) length(x) > 0))
		for(index in seq_along(s)) {
			tz = x@tracksCollection[[s[index]]]
			tz@tracks = tz@tracks[i[[s[index]]]]
			tz@tracksData = tz@tracksData[i[[s[index]]], ]
			# Write back the just processed Tracks element.
			x@tracksCollection[[s[index]]] = tz
		}
	} 
	# Drop data structure. Only relevant in case one single Tracks/Track element
	# have/has been selected. Multiple Tracks elements are always returned as
	# TracksCollection, independently of whether drop is true or false.
	if (drop && length(s) == 1) {
		if(is.list(i) && length(i[[s[1]]]) == 1)
			# Last [] is 1, since all but one Track elements have been sorted
			# out above.
			x@tracksCollection[[s]][1]
		else
			x@tracksCollection[[s]]
	}
	# Retain data structure, even if only one single Tracks/Track element
	# have/has been selected.
	else
		TracksCollection(x@tracksCollection[s], 
			x@tracksCollectionData[s,,drop=FALSE])
}

setMethod("[", "TracksCollection", subs.TracksCollection)

setMethod("[[", c("Track", "ANY", "missing"), 
	function(x, i, j, ...) {
		# TODO What if the attribute name coexists in both the data and
		# connections slot? Returning a list is inconvenient in the way that it
		# raises new design issues when making selections on objects of class
		# Tracks or TracksCollection: How to merge lists if tracks differ in
		# their "attribute state" (the passed attribute coexists in both slots,
		# exists in one slot only, does not exist)?
		if(i %in% names(x@data))
			x@data[[i]]
		else
			# Returns NULL if the attribute does not exist.
			x@connections[[i]]
	}
)

setMethod("[[", c("Tracks", "ANY", "missing"),
	function(x, i, j, ...) {
		do.call(c, lapply(x@tracks, function(t) t[[i]]))
	}
)

setMethod("[[", c("TracksCollection", "ANY", "missing"),
	function(x, i, j, ...) {
		do.call(c, lapply(x@tracksCollection, function(t) t[[i]]))
	}
)

setReplaceMethod("[[", c("Track", "ANY", "missing", "ANY"), 
	function(x, i, j, value) {
		if (i %in% names(x@connections)) {
			warning(paste("replacing", i, "in connections slot"))
			x@connections[[i]] = value
		} else
			x@data[[i]] = value
		x 	
	}
)

setReplaceMethod("[[", c("Tracks", "ANY", "missing", "ANY"), 
	function(x, i, j, value) {
		for(index in seq_along(x@tracks)) {
			if(i %in% names(x[index]@data)) {
				# "dim" (and with that "from" and "to") have to be reinitialized
				# each loop, since tracks might differ in their "attribute
				# state" (the passed attribute coexists in both slots, exists in
				# one slot only, does not exist).
				dim = sapply(x@tracks, function(t) dim(t)[["geometries"]])
				from = if(index == 1) 1 else cumsum(dim)[index-1] + 1
				to = cumsum(dim)[index]
				x@tracks[[index]]@data[[i]] = value[from:to]
			} else if(i %in% names(x[index]@connections)) {
				dim = sapply(x@tracks, function(t) dim(t)[["geometries"]]) - 1
				from = if(index == 1) 1 else cumsum(dim)[index-1] + 1
				to = cumsum(dim)[index]
				x@tracks[[index]]@connections[[i]] = value[from:to]
			}
		}
		x
	}
)

setReplaceMethod("[[", c("TracksCollection", "ANY", "missing", "ANY"), 
	function(x, i, j, value) {
		index = 1
		for(tz in seq_along(x@tracksCollection)) {
			for(t in seq_along(x[tz]@tracks)) {
				if(i %in% names(x[tz][t]@data)) {
					dim = do.call(c, lapply(x@tracksCollection, 
						function(tz) sapply(tz@tracks,
							function(t) dim(t)[["geometries"]])))
					from = if(index == 1) 1 else cumsum(dim)[index-1] + 1
					to = cumsum(dim)[index]
					x@tracksCollection[[tz]]@tracks[[t]]@data[[i]] = value[from:to]
				} else if(i %in% names(x[tz][t]@connections)) {
					dim = do.call(c, lapply(x@tracksCollection, 
						function(tz) sapply(tz@tracks,
							function(t) dim(t)[["geometries"]]))) - 1
					from = if(index == 1) 1 else cumsum(dim)[index-1] + 1
					to = cumsum(dim)[index]
					x@tracksCollection[[tz]]@tracks[[t]]@connections[[i]] = value[from:to]
				}
				index = index + 1
			}
		}
		x
	}
)

setMethod("$", "Track", function(x, name) x[[name]])

setMethod("$", "Tracks", function(x, name) x[[name]])

setMethod("$", "TracksCollection", function(x, name) x[[name]])

setReplaceMethod("$", "Track", 
	function(x, name, value) {
		x[[name]] = value
		x 
	}
)

setReplaceMethod("$", "Tracks", 
	function(x, name, value) {
		x[[name]] = value
		x 
	}
)

setReplaceMethod("$", "TracksCollection", 
	function(x, name, value) {
		x[[name]] = value
		x 
	}
)

# Provide stack, unstack and concatenate c() methods.

stack.TracksCollection = function (x, select, ...) {
	stopifnot(missing(select))
	splt = function(Tr) lapply(Tr@tracks, function(x) Tracks(list(x)))
	l = lapply(x@tracksCollection, function(x) splt(x))
	TracksCollection(do.call(c, l))
}

c.Tracks = function(...)
	Tracks(do.call(c, lapply(list(...), function(x) x@tracks)))

unstack.TracksCollection = function(x, form, ...) {
	TracksCollection(lapply(split(x@tracksCollection, form), 
		function(x) do.call(c, x)))
}

# approx coordinates and attributes for specific new time points, 
approxTrack = function(track, when, ..., n = 50, by, FUN = stats::approx, warn.if.outside = TRUE) {
	x = index(track)
	if (missing(when)) {
		if (missing(by))
			when = seq(min(x), max(x), length.out = n)
		else
			when = seq(min(x), max(x), by = by)
	} else if (warn.if.outside &&
		(min(when) < min(index(track)) || max(when) > max(index(track))))
			warning("approxTrack: approximating outside data range")
	p = apply(coordinates(track), 2, function(y) FUN(x, y, xout = when, n = n, ...)$y)
	d = data.frame(lapply(track@data, function(y) FUN(x, y, xout = when, n = n, ... )$y))
	Track(STIDF(SpatialPoints(p, CRS(proj4string(track))), when, d))
}
approxTracks = function(tr, ...) Tracks(lapply(tr@tracks, function(x) approxTrack(x,...)))
approxTracksCollection = function(tc, ...)
	TracksCollection(lapply(tc@tracksCollection, function(x) approxTracks(x,...)))
