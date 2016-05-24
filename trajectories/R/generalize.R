# Provide generalize methods.

generalize.Track <- function(t, FUN = mean, ..., timeInterval, distance, n, tol, toPoints) {
	if (sum(!c(missing(timeInterval), missing(distance), missing(n))) != 1)
		stop("exactly one parameter from (timeInterval, distance, n) has to be specified")
	if(!missing(timeInterval)) {
		origin = index(t@time)
		cut = cut(origin, timeInterval)
		segmentLengths = rle(as.numeric(cut))$lengths
	} 
	if (!missing(distance)) {
		# Total distances from each point to the first one.
		origin = c(0, cumsum(t@connections$distance))
		cut = floor(origin / distance)
		segmentLengths = rle(cut)$lengths
	} 
	if (!missing(n)) {
		dim = dim(t)["geometries"]
		if(n != 1 && dim / n > 1) {
			rep = floor((dim-n)/(n-1) + 1)
			mod = (dim-n) %% (n-1)
			if(mod == 0)
				segmentLengths = rep(n, rep)
			else
				segmentLengths = c(rep(n, rep), mod + 1)
		} else
			segmentLengths = dim
	} 
	# Update segment lengths to consider all segments for generalisation. In
	# case the cut-point falls between two points of the track to be
	# generalised, attach the next point to the current segment. If the cut-
	# point matches a point of the track, leave everything as is.
	toIndex = cumsum(segmentLengths)
	segmentLengths_ = integer()
	for(i in seq_along(segmentLengths)) {
		if (i == length(segmentLengths)
			|| (!missing(timeInterval) && origin[toIndex[i]] %in% seq(origin[1], origin[length(origin)], timeInterval))
			|| (!missing(distance) && origin[toIndex[i]] > 0 && origin[toIndex[i]] %% distance == 0)
			|| (!missing(n)))
			segmentLengths_[i] = segmentLengths[i]
		else { 
			segmentLengths_[i] = segmentLengths[i] + 1
			if(i == length(segmentLengths) - 1 && segmentLengths[i+1] == 1)
				break()
		}
	}
	segmentLengths = segmentLengths_
	# Aggregate over each segment.
	stidfs = list()
	endTime = NULL
	for(i in seq_along(segmentLengths)) {
		from = if(i == 1) 1 else tail(cumsum(segmentLengths[1:(i-1)]), n = 1) - (i-2)
		to = from + segmentLengths[i] - 1
		if(!missing(toPoints) && toPoints)
			sp = t@sp[(from+to)/2]
		else {
			l = Lines(list(Line(t@sp[from:to])), paste("L", i, sep = ""))
			sp = SpatialLines(list(l), proj4string = CRS(proj4string(t)))
			if(!missing(tol) && nrow(coordinates(sp)[[1]][[1]]) > 1) {
				if (!requireNamespace("rgeos", quietly = TRUE))
					stop("rgeos required for tolerance")
				sp = rgeos::gSimplify(spgeom = sp, tol = tol, 
					topologyPreserve = TRUE)
			}
		}
		time = t@time[from]
		if (is.null(endTime)) {
			endTime = t@endTime[to]
			tz = attr(endTime, "tzone")
		} else
			endTime = c(endTime, t@endTime[to])
		data = data.frame(lapply(t@data[from:to, , drop = FALSE], FUN, ...)) # EP added ...
		#stidfs = c(stidfs, STIDF(sp, time, data, t@endTime[to]))
		stidfs = c(stidfs, STIDF(sp, time, data))
	}
	stidf = do.call(rbind, stidfs)
	# Provide a workaround, since rbind'ing objects of class POSIXct as used
	# in the "endTime" slot of STIDF objects does not work properly.
	attr(endTime, "tzone") = tz
	stidf@endTime = endTime
	Track(stidf)
}

if(!isGeneric("generalize"))
	setGeneric("generalize", function(t, FUN = mean, ...)
		standardGeneric("generalize"))

setMethod("generalize", signature(t = "Track"), generalize.Track)

setMethod("generalize", signature(t = "Tracks"),
	function(t, FUN = mean, ...) {
		t@tracks = lapply(t@tracks,
			function(x) generalize(x, FUN, ...))
		t
	}
)

setMethod("generalize", signature(t = "TracksCollection"),
	function(t, FUN = mean, ...) {
		t@tracksCollection = lapply(t@tracksCollection,
			function(x) generalize(x, FUN, ...))
		t
	}
)
