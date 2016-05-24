create_new_data = function(obj) {
# Function that creates a new_data object if one is missing
	convex_hull = chull(coordinates(obj)[,1],coordinates(obj)[,2])
	convex_hull = c(convex_hull, convex_hull[1]) # Close the polygon
	d = Polygon(coordinates(obj)[convex_hull, ]) 
	new_data = spsample(d, 5000, type = "regular")
	gridded(new_data) = TRUE
	return(new_data)
}
