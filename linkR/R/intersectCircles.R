intersectCircles <- function(circle1, circle2){
	# http://math.stackexchange.com/questions/938701/how-to-find-the-intersection-points-of-two-circles-in-3d
	# TO ADD: ONE CIRCLE WITHIN THE OTHER, NO OVERLAP IN CIRCUMFERENCE

	# CURRENTLY ONLY WORKS WITH CIRCLES THAT ARE COPLANAR
	if(avectors(circle1$N, circle2$N)) stop(paste0("Currently only finds intersection of circles that are co-planar. Input circles are not coplanar; difference in angle between normal vectors: ", avectors(circle1$N, circle2$N)))

	# FIND DISTANCE BETWEEN CENTERS
	center_dist <- distPointToPoint(circle1$C, circle2$C)
	
	if(center_dist == 0){
		
		# CIRCLES PERFECTLY COINCIDENT
		if(circle1$R == circle2$R) return(list('type'='coincident'))

		# CIRCLES ARE CONCENTRIC
		if(circle1$R != circle2$R) return(list('type'='concentric'))
	}
	
	# CHECK THAT INTERSECTION IS POSSIBLE
	if(center_dist > circle1$R + circle2$R) return(list('type'='non-coincident'))

	# CHECK IF ONLY ONE SOLUTION EXISTS
	if(center_dist == circle1$R + circle2$R) return(list(circle1$C + uvector(circle2$C - circle1$C)*circle1$R, 'type'='one'))

	# FIND AREA OF TRIANGLE USING HERONS FORMULA
	p <- (circle1$R + circle2$R + center_dist) / 2
	area_sq <- p*(p - circle1$R)*(p - circle2$R)*(p - center_dist)

	# IF AREA IS NAN, MAYBE CIRCLE IS INSIDE THE OTHER (PLOT CONFIRMS ON ONE OCCASION)
	if(area_sq < 0) return(list('type'='inside'))	

	# FIND AREA
	area <- sqrt(area_sq)

	# FIND TRIANGLE HEIGHT
	h <- area / (0.5*center_dist)

	# FIND LENGTH FROM FIRST CENTER TO INTERSECTION MIDPOINT
	dist_c1m <- sqrt(circle1$R^2 - h^2)

	# PROJECT FROM FIRST CENTER TO INTERSECTION MIDPOINT
	inter_midpt <- circle1$C + uvector(circle2$C - circle1$C)*dist_c1m
	
	# DIRECTION VECTOR FROM INTERSECTION MIDPOINT TO INTERSECTIONS
	dir_vec <- uvector(circle2$C - circle1$C) %*% tMatrixEP(v=circle2$N, a=pi/2)

	# FIND TENTATIVE INTERSECTS
	intersects <- list(c(inter_midpt + dir_vec*h), c(inter_midpt + -dir_vec*h), 'type'='two')
	
	# MAKE SURE THAT INTERSECTION MIDPOINT WAS PROJECTED IN THE RIGHT DIRECTION
	if(abs(distPointToPoint(circle2$C, intersects[[1]]) - circle2$R) > 1e-8){
		
		# IF NOT, PROJECT IN OTHER DIRECTION
		inter_midpt <- circle1$C + -uvector(circle2$C - circle1$C)*dist_c1m
	
		# DIRECTION VECTOR FROM INTERSECTION MIDPOINT TO INTERSECTIONS
		dir_vec <- uvector(circle2$C - circle1$C) %*% tMatrixEP(v=circle2$N, a=pi/2)

		# FIND TENTATIVE INTERSECTS
		intersects <- list(c(inter_midpt + dir_vec*h), c(inter_midpt + -dir_vec*h), 'type'='two')
	}

	intersects
}