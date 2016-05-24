setClass("difftrack",
	slots=c(track1 ="Track", track2 = "Track", 
		conns1 = "SpatialLinesDataFrame", conns2 = "SpatialLinesDataFrame"),
)

## plots a difftrack
plot.difftrack <- function(x, y, ..., axes = TRUE) {  
	plot(x@track1@sp, col="red", ..., axes = axes)
	points(x@track2@sp, col="blue")
	lines(x@conns1)
	lines(x@conns2)
}

setMethod("plot", "difftrack", plot.difftrack)

## stcube for difftrack
stcube.difftrack <- function(x, showMap = FALSE, mapType = "osm", normalizeBy = "week", ..., y, z) {
  tracks <- Tracks(list(x@track1, x@track2))
  stcube(tracks, showMap = showMap, mapType = mapType, normalizeBy = normalizeBy, ...)  
  lines1 <- x@conns1@lines
  lines2 <- x@conns2@lines
  
  time1 <- x@conns1@data$time
  time1 <- normalize(time1, normalizeBy)
  time2 <- x@conns2@data$time
  time2 <- normalize(time2, normalizeBy)
  
  sapply(lines1, function(l) {
    coords <- coordinates(l)
    id <- as.numeric(l@ID)
    z <- time1[id]
    rgl::lines3d(coords[[1]][,1], coords[[1]][,2], z,  col = "red")    
  })
  sapply(lines2, function(l) {
    coords <- coordinates(l)
    id <- as.numeric(l@ID)
    z <- time2[id]
    rgl::lines3d(coords[[1]][,1], coords[[1]][,2], z, col = "blue")    
  })
}

setMethod("stcube", "difftrack", stcube.difftrack)
