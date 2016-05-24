closestCountry <- function(pt, proj = "+proj=longlat +datum=WGS84") {

	# check that provided proj4string is valid
	CRS(proj)

	if (proj != "+proj=longlat +datum=WGS84") {
		pt <- SpatialPoints(pt, proj4string=CRS(proj))
		pt <- spTransform(pt, CRS('+proj=longlat +datum=WGS84'))
	}
	if (is.matrix(pt) | is.data.frame(pt)) {
		pt <- as.numeric(pt)
	}
	d <- spDistsN1(worldPoints, pt, longlat = FALSE)
	w <- which.min(d)
	return(worldPointCountries[[w]])
}


