# Copyright 2004 Edzer Pebesma (copied from sp package)

"mt.point.in.polygon" <-
function(point.x, point.y, pol.x, pol.y) {
	.Call("R_point_in_polygon_mt", 
		as.numeric(point.x),
		as.numeric(point.y),
		as.numeric(pol.x),
		as.numeric(pol.y) 
		, PACKAGE = "maptools"
		)
}
