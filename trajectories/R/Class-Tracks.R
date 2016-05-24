# The first class is "Track", and contains a single track, or trip, followed by
# a person, animal or object. This means that consecutive location/time stamps
# are not interrupted by a period of substantially other activity. This object
# extends "STIDF", where locations and times as well as attributes measured at
# these locations and times, such as elevation, are stored. The function "Track"
# can now be used to build such an object from its components. The slot
# "connections" contains data about the segments between consecutive ST points.

setClass("Track",
	contains = "STIDF", # Locations, times and attribute data about the points.
	representation(connections = "data.frame"), 
	# With attribute data BETWEEN points: speed, direction etc.
	validity = function(object) {
		stopifnot(nrow(object@connections) + 1 == nrow(object@data))
		return(TRUE)
	}
)

directions_ll = function(cc, ll) {
	# cc a 2-column matrix with points, [x y] or [long lat]
	# ll a boolean, indicating longlat (TRUE) or not (FALSE)
	if (! ll) {
		dcc = matrix(apply(cc, 2, diff), ncol = ncol(cc))
		((atan2(dcc[,1], dcc[,2]) / pi * 180) + 360) %% 360
	} else {
		# longlat:
		# http://www.movable-type.co.uk/scripts/latlong.html
		# initial bearing:
		cc = cc * pi / 180
		lat1 = head(cc[,2], -1)
		lat2 = tail(cc[,2], -1)
		lon1 = head(cc[,1], -1)
		lon2 = tail(cc[,1], -1)
		dlon = lon2 - lon1
		az = atan2(sin(dlon)*cos(lat2),
			cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(dlon))
		((az / pi * 180) + 360) %% 360
	}
}

TrackStats = function(track) {
	duration = diff(as.numeric(index(track@time))) # seconds
	stopifnot(!any(duration == 0))
#	if (class(try(cc <- coordinates(track), silent=TRUE)) == "try-error" ||
#			!is.matrix(cc))
	if (!is(track@sp, "SpatialPoints"))
		data.frame(matrix(nrow = length(track@sp) - 1, ncol = 0)) # empty
	else {
		cc = coordinates(track@sp)
		ll = identical(is.projected(track), FALSE)
		distance = LineLength(cc, ll, FALSE) 
				# for sp 1.1-1: use spDists with segments = TRUE
		if (ll) # distance is in km, transform to m:
			distance = distance * 1000.0
		speed = distance / duration # per second
		direction = directions_ll(cc, ll)
		data.frame(distance = distance, duration = duration, 
			speed = speed, direction = direction)
	} 
}

# Computes segment lengths.

Track = function(track, df = fn(track), fn = TrackStats) {
	if (is(track, "STI") && !is(track, "STIDF"))
		track = STIDF(track@sp, track@time, data.frame(ones = rep(1, length(track))), track@endTime)
	duration = diff(as.numeric(index(track@time))) # seconds
	if (any(duration == 0)) {
		sel = (c(1, duration) != 0)
		n = sum(!sel)
		warning(paste("deselecting", n, "secondary point(s) of zero duration interval(s)"))
		if (sum(sel) < 2)
			stop("less than two unique time instances")
		track = Track(as(track, "STIDF")[sel,])
	}
	new("Track", track, connections = df)
}

# A collection of Track objects for single ID (car, person etc.).

setClass("Tracks", 
	representation(tracks = "list", tracksData = "data.frame"),
	validity = function(object) {
		stopifnot(all(sapply(object@tracks, function(x) is(x, "Track"))))
		stopifnot(nrow(object@tracksData) == length(object@tracks))
		stopifnot(length(object@tracks) > 0)
		stopifnot(!is.null(names(object@tracks)))
		stopifnot(identicalCRS(object@tracks))
		return(TRUE)
	}
)

TrackSummary = function(track) {
	ix = index(track@time)
	bb = bbox(track@sp)
	conn = track@connections
	data.frame(
		xmin = bb[1,1],
		xmax = bb[1,2],
		ymin = bb[2,1],
		ymax = bb[2,2],
		tmin = min(ix),
		tmax = max(ix),
		n = length(track@sp),
		distance = sum(conn$distance),
		medspeed = quantile(conn$speed, 0.5)
		# TODO Compute some mean direction?
	)
}

# Pre-computes elements of tracksData.

Tracks = function(tracks, 
		tracksData = data.frame(row.names=names(tracks)), fn = TrackSummary) {
	if (is.null(names(tracks)) && length(tracks) > 0)
		names(tracks) = paste("Track", 1:length(tracks), sep = "")
	new("Tracks", tracks = tracks, 
		tracksData = cbind(tracksData, do.call(rbind, lapply(tracks, fn))))
}

# Collection of Tracks for several IDs.
 
setClass("TracksCollection", 
	representation(tracksCollection = "list", 
		tracksCollectionData = "data.frame"),
	validity = function(object) {
		stopifnot(all(sapply(object@tracksCollection, class) == "Tracks"))
		stopifnot(length(object@tracksCollection) == 
			nrow(object@tracksCollectionData))
		stopifnot(length(object@tracksCollection) > 0)
		names = names(object@tracksCollection)
		stopifnot(!(is.null(names) || any(is.na(names))))
		stopifnot(identicalCRS(object@tracksCollection))
		return(TRUE)
	}
)

# unTracksCollection <- function(x, recursive=TRUE, use.names=TRUE) {
#   
#   TracksCollection(Tracks(lapply(x@tracksCollection,
#                                  function(tracksObj) {
#                                    unlist(tracksObj@tracks)
#                                  })))
#   
# }

TracksSummary = function(tracksCollection) {
	tc = tracksCollection
	df = data.frame(n = sapply(tc, function(x) length(x@tracks)))
	df$xmin = sapply(tc, function(x) min(x@tracksData$xmin))
	df$xmax = sapply(tc, function(x) max(x@tracksData$xmax))
	df$ymin = sapply(tc, function(x) min(x@tracksData$ymin))
	df$ymax = sapply(tc, function(x) max(x@tracksData$ymax))
	df$tmin = as.POSIXct(unlist(lapply(lapply(tc, function(x) x@tracksData$tmin), min)),
		origin = "1970-01-01", tz=attr(tc[[1]]@tracks[[1]]@time, "tz"))
		# do.call(c, lapply(lapply(tc, function(x) x@tracksData$tmin), min)) # reported by RH
	df$tmax = as.POSIXct(unlist(lapply(lapply(tc, function(x) x@tracksData$tmax), max)),
		origin = "1970-01-01", tz=attr(tc[[1]]@tracks[[1]]@time, "tz"))
		# do.call(c, lapply(lapply(tc, function(x) x@tracksData$tmax), max)) # reported by RH
	row.names(df) = names(tracksCollection)
	df
}

TracksCollection = function(tracksCollection, tracksCollectionData = NULL,
	fn = TracksSummary) {
	if (is.null(names(tracksCollection)))
		names(tracksCollection) = paste("Tracks", 1:length(tracksCollection), 
		sep = "")
	ts = TracksSummary(tracksCollection)
	if (is.null(tracksCollectionData))
		tracksCollectionData = ts
	else
		tracksCollectionData = cbind(tracksCollectionData, ts)
	new("TracksCollection", tracksCollection = tracksCollection, 
		tracksCollectionData = tracksCollectionData)
}
