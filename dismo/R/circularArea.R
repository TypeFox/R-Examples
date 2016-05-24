
.circularArea <- function(p, radius, lonlat=FALSE) {
	requireNamespace('rgeos')
	x <- circles(p, d=radius, lonlat=lonlat)
	if (lonlat) {
	
	} else {
		rgeos::gArea(x@polygons) / 1000000
	}
}
