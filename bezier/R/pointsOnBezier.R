pointsOnBezier <- function(p, n = NULL, method = 'evenly_spaced', t1 = 0, t2 = NULL, deg = NULL, max.dist = NULL, max.dist.factor = 0.1, relative.min.slope = 1e-7, absolute.min.slope = 0, sub.relative.min.slope = 1e-4, sub.absolute.min.slope = 0, print.progress = FALSE){

	# IF N IS NULL AND MAX.DIST IS DEFINED, METHOD IS MAX_DIST
	if(is.null(n) && !is.null(max.dist)) method <- 'max_dist'

	# CHECK THAT N AND MAX.DIST ARE NOT BOTH NULL
	if(tolower(method) %in% c('max_dist', 'evenly_spaced') && is.null(n) && is.null(max.dist)) stop(paste0("With the ", method, " method, n and max.dist cannot both be NULL. Specify either the number of points to generate (n) or the maximum distance between consecutive points (max.dist)."))

	# IF P IS A LIST WITH PARAMETERS AS SEPARATE DIMENSIONS, CONVERT TO MATRIX, ELEVATING LOWER DEGREES, IF NECESSARY
	if(is.list(p)){
	
		# FIND MAXIMUM DEGREE
		max_deg <- max(unlist(lapply(p, 'length'))) - 1
		
		# ELEVATE EACH PARAMETRIC BEZIER TO MAXIMUM DEGREE
		for(i in 1:length(p)) p[[i]] <- elevateBezierDegree(p[[i]], deg = max_deg)
		
		# UNLIST AND CONVERT TO MATRIX
		p <- matrix(unlist(p), ncol=length(p))
	}

	# IF DEGREE IS NULL, DEFAULT IS NUMBER OF ROWS IN P MINUS ONE
	if(is.null(deg)) if(is.vector(p)){deg <- length(p) - 1}else{deg <- nrow(p) - 1}

	# IF T2 IS NULL, FIND TOTAL T BASED ON BEZIER DEGREE
	if(is.null(t2)) if(is.vector(p)){t2 <- (length(p) - 1) / deg}else{t2 <- (nrow(p) - 1) / deg}

	# SPLIT BY METHOD
	if(method == 'adjoining'){

		# GET BEZIER ARC LENGTH
		bezier_arc_length <- bezierArcLength(p, t1=t1, t2=t2, deg=deg, relative.min.slope=1)$arc.length

		# SET INTIAL S1
		s1 <- t1

		# SET INITIAL INTERVAL
		iter <- (t2-t1) / bezier_arc_length
		
		# SET INTIAL T AND ERROR VECTORS
		points_on_bezier <- matrix(p[1, ], nrow=1, ncol=ncol(p))
		error <- rep(NA, 0)
		t2_vector <- c(0)
		
		while(s1 < t2){
			
			# GET POINT POSITION
			point_position <- bezier(t=s1 + iter, p=p, deg=deg)[1, ]
			point_position_round <- round(point_position)
			
			# GET DISTANCE FROM PREVIOUS POINT
			dist_to_prev <- abs(points_on_bezier[nrow(points_on_bezier), ] - point_position_round)
			
			# IF DISTANCE IS TOO SMALL, INCREASE INTERVAL
			if(sum(dist_to_prev) == 0){
				iter <- iter + iter
				next
			}

			# IF DISTANCE IS TOO LARGE, DECREASE INTERVAL
			if(dist_to_prev[1] > 1 || dist_to_prev[2] > 1){
				iter <- iter - iter*0.4
				next
			}
			
			# IF POINTS ARE ADJOINING SAVE POINT AND DISTANCE FROM ROUNDED TO ACTUAL POSITION
			points_on_bezier <- rbind(points_on_bezier, point_position_round)			
			error <- c(error, sqrt(sum((point_position - point_position_round)^2)))
			
			# SAVE PARAMETRIC VALUE
			if(s1 + iter > t2){
				t2_vector <- c(t2_vector, t2)
			}else{
				t2_vector <- c(t2_vector, s1 + iter)
			}

			# UPDATE STEP VARIABLE
			s1 <- s1 + iter
		}

		rownames(points_on_bezier) <- NULL
		return(list(points = points_on_bezier, error = error, t = t2_vector))
	}

	# GET BEZIER ARC LENGTH
	bezier_arc_length <- bezierArcLength(p, t1=t1, t2=t2, deg=deg, relative.min.slope=relative.min.slope, absolute.min.slope=absolute.min.slope)$arc.length

	# RAPID FIND POINTS ON BEZIER BY CONSTRAINING THE MAXIMUM DISTANCE BETWEEN ADJACENT POINTS
	if(method == 'max_dist'){
		
		# SET INTIAL S1
		s1 <- t1
		
		# CHECK THAT MAXIMUM DISTANCE DOES NOT EXCEED HALF ARC LENGTH
		if(max.dist > bezier_arc_length/2) stop(paste0("Specified maximum distance (", max.dist, ") exceeds half the total arc length (", round(bezier_arc_length, 4), ")"))
		
		# SET INTIAL INTERVAL
		iter <- (t2-t1)/(bezier_arc_length/max.dist)

		# SET INTIAL T AND ERROR VECTORS
		t2_vector <- rep(0, 1)
		error <- rep(NA, 0)

		while(s1 < t2){

			# IF INTERVAL EXCEEDS T2, REPLACE WITH INTERVAL BETWEEN S1 AND T2
			if(s1 + iter > t2) iter <- t2 - s1

			# IF DISTANCE IS LESS THAN INPUT MAX DIST, INCREASE ITER UNTIL DISTANCE JUST EXCEEDS MAX DIST AND TAKE PENULTIMATE ITER
			if(iter != t2 - s1 && sqrt(sum((bezier(t=s1, p, deg=deg) - bezier(t=s1 + iter, p, deg=deg))^2)) < max.dist){
				while(sqrt(sum((bezier(t=s1, p, deg=deg) - bezier(t=s1 + iter, p, deg=deg))^2)) < max.dist) iter <- iter + iter*max.dist.factor
				iter <- iter - iter*max.dist.factor
			}

			# IF DISTANCE IS GREATER THAN INPUT MAX DIST, REDUCE INTERVAL UNTIL DISTANCE IS LESS THAN MAX DIST
			while(sqrt(sum((bezier(t=s1, p, deg=deg) - bezier(t=s1 + iter, p, deg=deg))^2)) > max.dist) iter <- iter - iter*max.dist.factor

			# SAVE VALUES TO VECTORS
			t2_vector <- c(t2_vector, s1 + iter)
			error <- c(error, max.dist - sqrt(sum((bezier(t=s1, p, deg=deg) - bezier(t=s1 + iter, p, deg=deg))^2)))

			# NEXT S1 IS S1 PLUS CURRENT ITER
			s1 <- s1 + iter
		}

		return(list(points = bezier(t=t2_vector, p=p, deg=deg), error = error, t = t2_vector))
	}

	## EVEN-SPACING METHOD

	# GET LENGTH OF EACH SEGMENT
	seg_length <- bezier_arc_length/(n-1)

	# SET INTIAL T AND ERROR VECTORS
	t2_vector <- c(rep(t1, n-1), t2)
	error <- rep(0, n)
	
	if(print.progress) cat('\nn of ', n, ': 1 ', sep='')

	for(i in 2:(n-1)){
		if(print.progress) cat(i, ' ', sep='')

		# INITIAL GUESS FOR PARAMETER TO OPTIMIZE
		par <- t2_vector[i-1] + (t2 - t2_vector[i-1])*(1/(n - i + 1))

		# FIND T FOR WHICH BEZIER LENGTH IS EQUAL TO LENGTH FROM END OF SEGMENT TO START
		optim_r <- optim(par = par, compareBezierArcLength, p = p, l = seg_length*(i-1), t1 = t1, deg = deg, relative.min.slope = sub.relative.min.slope, absolute.min.slope = sub.absolute.min.slope, method = "Brent", lower = t2_vector[i-1], upper = t2)

		# SAVE OPTIMUM VALUE OF T2 FOR WHICH ARC LENGTH EQUALS SPECIFIED SEGMENT LENGTH
		t2_vector[i] <- optim_r$par

		# DIFFERENCE BETWEEN MINIMIZED ARC LENGTH AND EXPECTED SEGMENT LENGTH
		error[i] <- optim_r$value
	}

	if(print.progress) cat('', n, '\n', sep='')

	# RETURN POINTS AND ERROR VECTORS
	list(points = bezier(t=t2_vector, p=p, deg=deg), error = error, t = t2_vector)
}