findEpipolarTangencyAngles <- function(curve_points, cal.coeff, window.tangency){

	avectors_save <- rep(NA, nrow(curve_points))

	for(i in 2:(nrow(curve_points)-1)){

		# FIND SELF EPIPOLAR LINE
		self_epipolar_slope <- abs(dltEpipolarLine(p=curve_points[i, ], cal.coeff1=cal.coeff, self=TRUE)$m)

		# VERIFY THAT REFERENCE POINT IS ON SELF-EPIPOLAR (ZERO DISTANCE)
		#sel <- dltEpipolarLine(p=curve_points[i, ], cal.coeff1=cal.coeff, self=TRUE)
		#dptl <- distancePointToLine(p=curve_points[i, ], sel$l1, sel$l2)

		# GET INDICES FOR DETERMINING SLOPE
		s_min <- max(1, i - window.tangency)
		s_max <- min(i + window.tangency, nrow(curve_points))

		# FIND SLOPE COMPONENTS
		slope <- rowSums(cbind(curve_points[s_min, ] - curve_points[i, ], curve_points[i, ] - curve_points[s_max, ]))

		# FIND CURVE SLOPE
		if(slope[1] == 0){curve_slope <- 10^4}else{curve_slope <- abs(slope[2]/slope[1])}
		
		# FIND ANGLE BETWEEN SLOPES IN RADIANS
		a_vectors <- avectors(c(1, curve_slope), c(1, self_epipolar_slope))
		avectors_save[i] <- a_vectors
	}

	avectors_save
}
