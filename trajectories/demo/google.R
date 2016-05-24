# read a file obtained by "Export to kml" from https://maps.google.com/locationhistory/
filename = system.file("history-06-09-2015.kml", package="trajectories")

library(XML)
kml <- xmlToList(filename)

tr = kml$Document$Placemark$Track
cc = which(names(tr) == "coord")
coord = t(sapply(kml$Document$Placemark$Track[cc], function(x) scan(text = x, quiet = TRUE)))[,1:2]
when = which(names(tr) == "when")
# convert the "-07:00" into " -0700" with sub:
#time = strptime(sub("-08:00$", " -0800", unlist(kml$Document$Placemark$Track[when])),
time = strptime(sub("([+\\-])(\\d\\d):(\\d\\d)$", "\\1\\2\\3", 
	unlist(kml$Document$Placemark$Track[when])), "%Y-%m-%dT%H:%M:%OS %z")

library(sp)
library(spacetime)
library(trajectories)
track = Track(STI(SpatialPoints(coord, CRS("+proj=longlat +ellps=WGS84")), time))
summary(track)
head(as(track, "data.frame"))
plot(track, axes = TRUE)

# the following will try to read your complete Locationhistory JSON dump:

library(jsonlite)
system.time(x <- fromJSON("Location History/LocationHistory.json")$locations)
object.size(x)
sapply(x, class)
x$time = as.POSIXct(as.numeric(x$timestampMs)/1000, origin = "1970-01-01")

x$lat = loc$latitudeE7 / 1e7
x$lon = loc$longitudeE7 / 1e7

a = x$activitys
types = unique(unlist(lapply(a, function(x) if (is.null(x[[2]])) "null" else x[[2]][[1]]$type)))
types = types[types != "null"] # NULL entries
getType = function(x, Type) {
	if (is.null(x)) 
		as.numeric(NA) 
	else {
		s = subset(x$activities[[1]], type == Type)$confidence
		if (length(s) == 0) 
			s = 0.0
		s[1]  # there might be more than one, we're now ignoring all others!
	}
}

for (Tp in types) # untangle:
	x[[Tp]] = sapply(a, getType, Type = Tp)

a = apply(x[,types], 1, function(x) if(all(is.na(x))) NA else which.max(x))
x$activity = factor(types[a])
x$confidence = apply(x[,types], 1, function(x) if(all(is.na(x))) NA else x[which.max(x)])
x$activitys = NULL # remove the complex list

library(sp)
loc.sp = x
coordinates(loc.sp) = ~lon+lat
proj4string(loc.sp) = CRS("+proj=longlat +datum=WGS84")

library(spacetime)
library(trajectories)
tr = Track(STIDF(geometry(loc.sp), loc.sp$time, loc.sp@data))
