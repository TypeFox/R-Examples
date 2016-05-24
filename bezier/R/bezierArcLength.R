bezierArcLength <- function(p, t1 = 0, t2 = NULL, deg=NULL, relative.min.slope = 1e-6, absolute.min.slope = 0, max.iter = 20, n = NULL){
	#cat('call bezierArcLength\n')

	# IF T1 IS EQUAL TO T2 RETURN ZERO LENGTH
	if(!is.null(t2) && t1 == t2){
		l <- list(arc.length=0, slope.break=NA, n=NA, break.cause=NA, n.iter=NA)
		class(l) <- 'bezierArcLength'
		return(l)
	}

	# IF DEGREE IS NULL, DEFAULT IS NUMBER OF ROWS IN P MINUS ONE
	if(is.null(deg)) if(is.vector(p)){deg <- length(p) - 1}else{deg <- nrow(p) - 1}

	# IF T2 IS NULL, FIND TOTAL T BASED ON BEZIER DEGREE
	if(is.null(t2)) if(is.vector(p)){t2 <- (length(p) - 1) / deg}else{t2 <- (nrow(p) - 1) / deg}

	# IF T2 EXCEEDS 1 (INPUT IS BEZIER SPLINE), LOOP THROUGH EACH SUB-SECTION OF SPLINE
	if(t2 - t1 > 1){

		# FIND NUMBER OF SUB-SECTIONS
		seg <- ceiling(t2) - floor(t1)
		
		# SET INITIAL T1 AND T2 SUB INTERVAL
		t1_sub <- t1
		
		# RETURN VECTORS AND VALUES
		arc_length <- 0
		slope_break <- rep(NA, seg)
		n_return <- rep(NA, seg)
		break_cause <- rep(NA, seg)
		n_iter <- rep(NA, seg)

		for(i in 1:seg){

			# FIND T2 OF SUB-INTERVAL 
			t2_sub <- min(t2, floor(t1_sub + 1))
			
			# CALL BEZIER ARC LENGTH FOR EACH SUB-INTERVAL
			bezier_arc_length <- bezierArcLength(p, t1 = t1_sub, t2 = t2_sub, deg=deg, relative.min.slope = relative.min.slope, absolute.min.slope = absolute.min.slope, n = n)

			# PUT RETURN VALUES INTO VECTORS
			arc_length <- arc_length + bezier_arc_length$arc.length
			slope_break[i] <- bezier_arc_length$slope.break
			n_return[i] <- bezier_arc_length$n
			break_cause[i] <- bezier_arc_length$break.cause
			n_iter[i] <- bezier_arc_length$n.iter

			# SET CURRENT T2 TO T1 OF NEXT SUB-INTERVAL
			t1_sub <- t2_sub
		}
		
		l <- list(arc.length=arc_length, slope.break=slope_break, n=n_return, break.cause=break_cause, n.iter=n_iter)
		class(l) <- 'bezierArcLength'
		return(l)
	}

	# IF NUMBER OF POINTS IS SPECIFIED, FIND SUM OF INTERPOINT DISTANCES AND RETURN
	if(!is.null(n)){
		
		# GET POINTS
		b <- bezier(t=seq(t1, t2, length=n), p=p, deg=deg)
		
		# FIND INTERPOINT DISTANCES
		if(is.vector(b)){d <- abs(b[2:length(b)] - b[1:(length(b)-1)])}else{d <- sqrt(rowSums((b[2:nrow(b), ] - b[1:(nrow(b)-1), ])^2))}

		# SUM INTERPOINT DISTANCES
		l <- list(arc.length=sum(d), n=n, slope.break=NA, n.iter=NA, break.cause=NA)
		class(l) <- 'bezierArcLength'
		return(l)
	}

	# NUMBER OF POINTS TO USE FOR SLOPE CALCULATION (USING FEWER SLOWS FUNCTION ~ 2X, WORSE PREDICTIONS?)
	slope_point_num <- 10

	# STEPS FOR INITIAL LENGTH CALCULATION -- FOR MOSTLY STRAIGHT CURVES WITH SHARP CONVOLUTION, NEED TO START WITH GREATER NUMBER
	t_lengths <- round(seq(5, 15, length=slope_point_num))

	# MATRIX FOR SAVING CUMULATIVE BEZIER LENGTHS
	b_len <- matrix(NA, nrow=0, ncol=2)

	# GIVE VALUE TO OUTPUT SLOPE ACHIEVED
	slope.break <- NA

	# BREAK FROM MAXIMUM ITERATIONS IF NOTHING ELSE
	break.cause <- paste('maximum number of iterations (', max.iter,') exceeded', sep='')
	
	for(i in 1:max.iter){

		# GET SUM OF BEZIER INTERPOINT DISTANCES FOR SPECIFIED LENGTHS
		for(j in 1:length(t_lengths)){

			# GET POINTS
			b <- bezier(t=seq(t1, t2, length=t_lengths[j]), p=p, deg=deg)
		
			# FIND INTERPOINT DISTANCES
			if(is.vector(b)){d <- abs(b[2:length(b)] - b[1:(length(b)-1)])}else{d <- sqrt(rowSums((b[2:nrow(b), ] - b[1:(nrow(b)-1), ])^2))}

			# SUM INTERPOINT DISTANCES AND ADD TO MATRIX
			b_len <- rbind(b_len, cbind(t_lengths[j], sum(d)))
		}

		# IF LENGTH IS UNIFORM FOR INITIAL LENGTH SET, CHECK THAT FUNCTION DOES NOT MISS SOME LATER CONVOLUTION
		if(i == 1 && sd(b_len[, 2]) == 0){

			# SET VECTOR OF LENGTHS TO TEST FOR UNIFORMITY
			t_lengths_proj <- 200

			# GET SUM OF BEZIER INTERPOINT DISTANCES FOR SPECIFIED LENGTHS
			for(j in 1:length(t_lengths_proj)){

				# GET POINTS
				b <- bezier(t=seq(t1, t2, length=t_lengths_proj[j]), p=p, deg=deg)

				# FIND INTERPOINT DISTANCES
				if(is.vector(b)){d <- abs(b[2:length(b)] - b[1:(length(b)-1)])}else{d <- sqrt(rowSums((b[2:nrow(b), ] - b[1:(nrow(b)-1), ])^2))}

				# SUM INTERPOINT DISTANCES AND ADD TO MATRIX
				b_len <- rbind(b_len, cbind(t_lengths_proj[j], sum(d)))
			}

			# RETURN IF STANDARD DEVIATION WITH NEW POINTS IS STILL ZERO
			if(sd(b_len[, 2]) == 0){
				l <- list(arc.length=b_len[nrow(b_len), 2], slope.break=0, n=b_len[nrow(b_len), 1], break.cause='uniform arc length for initial set', n.iter=NA)
				class(l) <- 'bezierArcLength'
				return(l)
			}
		}

		# FIND SLOPE OF LAST TEN POINTS
		slope <- summary(lm(b_len[(nrow(b_len)-slope_point_num):nrow(b_len), 2] ~ b_len[(nrow(b_len)-slope_point_num):nrow(b_len), 1]))$coefficients[2, 1]

		# BREAK IF SLOPE RELATIVE TO CURRENT CURVE LENGTH IS LESS THAN SPECIFIED MINIMUM
		if(slope/max(b_len[, 2]) < relative.min.slope){slope.break <- slope/max(b_len[, 2]);break.cause <- 'relative minimum slope reached';break}

		# BREAK IF ABSOLUTE SLOPE IS LESS THAN SPECIFIED MINIMUM
		if(slope < absolute.min.slope){slope.break <- slope;break.cause <- 'absolute minimum slope reached';break}

		# SET MAGNITUDE OF DIFFERENCE IN SLOPES TO ZERO
		mag_diff_slope <- c(0, 0)

		# ORDERS OF MAGNITUDE CURRENT SLOPE IS GREATER THAN SPECIFIED TOLERANCES
		if(absolute.min.slope > 0) mag_diff_slope[1] <- log(abs(slope)/absolute.min.slope, 10)
		if(relative.min.slope > 0) mag_diff_slope[2] <- log((abs(slope)/b_len[nrow(b_len), 2])/relative.min.slope, 10)

		# GET MAXIMUM
		mag_diff_slope <- min(mag_diff_slope)
		
		# FIND START OF NEXT SET OF LENGTHS
		next_n <- round(max(t_lengths) + max(t_lengths)*mag_diff_slope)
	
		# GET NEW LENGTHS
		t_lengths <- seq(next_n, next_n + slope_point_num)

		#cat('slope: ', slope, '\n')
		#cat('length: ', b_len[nrow(b_len), 2], '\n')
		#cat('\n')
	}
	
	if(break.cause == paste('maximum number of iterations (', max.iter,') exceeded', sep='')) slope.break <- slope

	#cat('slope.break: ', slope.break, '\n')
	#print(b_len)

	l <- list(arc.length=max(b_len[, 2]), slope.break=slope.break, n=b_len[which.max(b_len[, 2]), 1], break.cause=break.cause, n.iter=i)
	class(l) <- 'bezierArcLength'
	l
}