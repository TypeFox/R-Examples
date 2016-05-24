setCoordinateAxes <- function(pt1, pt2=NULL, pt3=NULL){

	# Returns 3 unit vectors (xu, yu, zu), each perpendicular to the others
	# xu is the uvector from pt1 to pt3 (the only vector which will always parallel to a line defined by two of the original points)
	# zu is a uvector perpendicular to the plane defined by pt1, pt2 and pt3
	# yu is a uvector perpendicular to zu and xu

	if(is.null(pt2)){
		pt_m <- pt1
		pt1 <- pt_m[1, ]
		pt2 <- pt_m[2, ]
		pt3 <- pt_m[3, ]
	}
	
	v2 <- pt3 - pt1
	zu <- uvector(cprod(pt2 - pt1, v2))
	xu <- uvector(v2)
	yu <- cprod(zu, xu)
	matrix(c(xu, yu, zu), 3, 3, byrow=T) # x,y,z vectors by row
}