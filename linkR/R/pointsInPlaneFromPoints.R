pointsInPlaneFromPoints <- function(p, d, n, compare=NULL){
	
	#print(list(p=p, d=d, n=n, compare=compare))

	# FIND ORTHOGONAL PROJECTION OF POINT1 ONTO PLANE
	center1 <- pointPlaneProj(p[1, ], p[2, ], n)

	# CIRCLE IN PLANE AT NORMAL FROM FIRST POINT
	circle1 <- defineCircle(center=center1, nvector=n, radius=sqrt(d[1]^2 - distPointToPoint(p[1, ], center1)^2), redefine_center=FALSE)

	# FIND ORTHOGONAL PROJECTION OF POINT ONTO PLANE
	center2 <- pointPlaneProj(p[5, ], p[4, ], n)

	# CIRCLE IN PLANE AT NORMAL FROM FIRST POINT
	circle2 <- defineCircle(center=center2, nvector=n, radius=sqrt(d[4]^2 - distPointToPoint(p[5, ], center2)^2), redefine_center=FALSE)

	# MOVE CIRCLE1 TO CIRCLE2 PLANE
	center12 <- pointPlaneProj(q=center1, p=center2, n=n)
	center21 <- pointPlaneProj(q=center2, p=center1, n=n)

	# FIND INTERSECTION POINTS OF CIRCLE1 ON CIRCLE2
	circle1$C <- center12
	iC12 <- intersectCircles(circle1, circle2)

	# FIND INTERSECTION POINTS OF CIRCLE2 ON CIRCLE1
	circle1$C <- center1
	circle2$C <- center21
	iC21 <- intersectCircles(circle1, circle2)
	
	# IF CIRCLES ARE PERFECTLY OVERLAPPING, RETURN CURRENT POINTS
	if(iC12[['type']] == 'coincident' || iC21[['type']] == 'coincident') return(p[2:4, ])

	# IF ONE CIRCLE IS INSIDE THE OTHER, NO OVERLAP, RETURN NA
	if(iC12[['type']] == 'inside' || iC21[['type']] == 'inside') return(NULL)

	# IF MORE THAN ONE POINT IS RETURNED AND COMPARE POINT PROVIDED, FIND CLOSEST POINT TO COMPARE
	min_idx12 <- min_idx21 <- 1
	if(length(iC12) > 1 && !is.null(compare)){
		min_idx12 <- which.min(c(distPointToPoint(compare[4, ], iC12[[1]]), distPointToPoint(compare[4, ], iC12[[2]])))
	}
	if(length(iC21) > 1 && !is.null(compare)){
		min_idx21 <- which.min(c(distPointToPoint(compare[2, ], iC21[[1]]), distPointToPoint(compare[2, ], iC21[[2]])))
	}

	return(rbind(
		iC21[[min_idx21]],
		pointPlaneProj(q=iC21[[min_idx21]], p=p[3, ], n=n),
		iC12[[min_idx12]]
	))
}