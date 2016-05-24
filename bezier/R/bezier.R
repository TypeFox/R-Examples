bezier <- function(t, p, start = NULL, end = NULL, deg=NULL){

	# IF P IS A LIST WITH PARAMETERS AS SEPARATE DIMENSIONS, CONVERT TO MATRIX, ELEVATING LOWER DEGREES, IF NECESSARY
	if(is.list(p)){
	
		# FIND MAXIMUM DEGREE
		max_deg <- max(unlist(lapply(p, 'length'))) - 1
		
		# ELEVATE EACH PARAMETRIC BEZIER TO MAXIMUM DEGREE
		for(i in 1:length(p)) p[[i]] <- elevateBezierDegree(p[[i]], deg = max_deg)
		
		# IF INPUT IS VECTOR RETURN AS VECTOR
		if(length(p) > 1){return_vector <- FALSE}else{return_vector <- TRUE}

		# UNLIST AND CONVERT TO MATRIX
		p <- matrix(unlist(p), ncol=length(p))
	}else{
		# IF INPUT IS VECTOR RETURN AS VECTOR
		if(is.matrix(p)){return_vector <- FALSE}else{return_vector <- TRUE}
	}

	# IF START POINT IS PROVIDED, APPEND TO START OF P VECTOR/MATRIX
	if(!is.null(start)) if(length(start) == 1){p <- c(start, p)}else{p <- rbind(start, p)}		

	# IF END POINT IS PROVIDED, APPEND TO END OF P VECTOR/MATRIX
	if(!is.null(end)) if(length(end) == 1){p <- c(p, end)}else{p <- rbind(p, end)}

	# CONVERT INPUT PARAMETER POINTS TO MATRIX
	if(!is.matrix(p)) p <- matrix(p)

	# IF DEGREE IS NULL, DEFAULT IS NUMBER OF ROWS IN P MINUS ONE
	if(is.null(deg)) deg <- nrow(p) - 1

	# USE DEGREE TO FIND NUMBER OF CONCATENATED SETS OF BEZIER CURVES
	num_s <- (nrow(p) - 1) / deg

	# CHECK THAT NUMBER OF PARAMETER ROWS MATCHES DEGREE
	if(num_s - floor(num_s) > 0) stop("Number of rows in parameter matrix do not match input degree.")

	# MATRIX FOR BEZIER PARAMETRIC POINTS
	b <- matrix(0, nrow=length(t), ncol=ncol(p))

	# SET ROW INDICES FOR EACH BEZIER SPLINE SEGMENT
	seg <- matrix(1, nrow=num_s, ncol=2)
	for(i in 1:num_s) seg[i, ] <- c((i-1)*deg + 1, (i-1)*deg + 1 + deg)
	
	for(j in 1:length(t)){
		
		# GET SPLINE SEGMENT NUMBER
		if(t[j] == 0){s <- 1}else{s <- ceiling(t[j])}

		# IF SEGMENT NUMBER EXCEEDS NUMBER OF SEGMENTS SET TO LAST (FOR PARAMETRIC VALUES JUST OVER THE CUT-OFF)
		if(s > num_s) s <- num_s

		# GET CORRESPONDING ROWS OF P MATRIX
		p_sub <- matrix(p[seg[s, 1]:seg[s, 2], ], nrow=deg+1, ncol=ncol(p))

		# FIND BEZIER POINTS FROM GENERALIZED BEZIER FORMULA
		b[j, ] <- colSums(choose(deg, 0:deg)*((1 - (t[j] - s + 1))^(deg - 0:deg))*(t[j] - s + 1)^(0:deg)*p_sub)
	}

	# IF INPUT WAS VECTOR RETURN AS VECTOR
	if(return_vector) return(as.vector(b))

	b
}