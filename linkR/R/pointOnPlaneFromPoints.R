pointOnPlaneFromPoints <- function(p, n, p1, d1, p2, d2, compare=NULL){
	
	#print(list(p=p, n=n, p1=p1, d1=d1, p2=p2, d2=d2))

	# FIND ORTHOGONAL PROJECTION OF POINT ONTO PLANE
	center1 <- pointPlaneProj(p1, p, n)

	# CIRCLE IN PLANE AT NORMAL FROM FIRST POINT
	circle1 <- defineCircle(center=center1, nvector=n, radius=sqrt(d1^2 - distPointToPoint(p1, center1)^2), redefine_center=FALSE)

	# FIND ORTHOGONAL PROJECTION OF POINT ONTO PLANE
	center2 <- pointPlaneProj(p2, p, n)

	# CIRCLE IN PLANE AT NORMAL FROM FIRST POINT
	circle2 <- defineCircle(center=center2, nvector=n, radius=sqrt(d2^2 - distPointToPoint(p2, center2)^2), redefine_center=FALSE)

	# FIND INTERSECTION POINTS OF CIRCLES
	iC <- intersectCircles(circle1, circle2)
	
	# IF CIRCLES ARE PERFECTLY OVERLAPPING, RETURN CURRENT POINT
	if(iC[['type']] == 'coincident') return(list(p))

	# IF ONE CIRCLE IS INSIDE THE OTHER, NO OVERLAP, RETURN NA
	if(iC[['type']] == 'inside') return(list(NA))

	# IF MORE THAN ONE POINT IS RETURNED AND COMPARE POINT PROVIDED, FIND CLOSEST POINT TO COMPARE
	if(length(iC) > 1 && !is.null(compare)){
		
		min_idx <- which.min(c(distPointToPoint(compare, iC[[1]]), distPointToPoint(compare, iC[[2]])))
		
		return(list(iC[[min_idx]]))
	}
	
	return(iC)
}