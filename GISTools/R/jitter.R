jitter.points <- function(pts,scl) {
	x = coordinates(pts)
	x =  x + rnorm(length(x),0,scl) 
	res = SpatialPoints(x)
	proj4string(res)=CRS(proj4string(pts))
	if (class(pts)=="SpatialPointsDataFrame") {
		res = SpatialPointsDataFrame(res,data.frame(pts))}
	return(res) }

bstrap.points <- function(pts) {
	x = coordinates(pts)
	x = x[sample(nrow(x),replace=TRUE),]
	res = SpatialPoints(x)
	proj4string(res)=CRS(proj4string(pts))
	if (class(pts)=="SpatialPointsDataFrame") {
		res = SpatialPointsDataFrame(res,data.frame(pts))}
	return(res) }
	
poly.labels <- function(polys) {
	pts = SpatialPoints(t(sapply(slot(polys,'polygons'), function(x) slot(x,'labpt'))))
	proj4string(pts) = CRS(proj4string(polys))
	pts }

