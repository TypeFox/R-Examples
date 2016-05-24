matchCurvePoints <- function(curve_ref, curve_1, cal_coeff, min.direct.tangency = 25, 
	min.fill.tangency = 0, epi.err.weight = 1, rec.err.weight = 0, window.tangency = 5, 
	direct.min = 4){

	# MATRIX OF POINTS HOMOLOGOUS TO REF CURVE
	curve_1_hom <- matrix(NA, nrow=nrow(curve_ref), ncol=ncol(curve_ref))

	# FIND EPIPOLAR TANGENCY ANGLES FOR BOTH VIEWS
	tan_angles <- findEpipolarTangencyAngles(curve_ref, cal_coeff[, 1:2], window.tangency)

	# MATRIX OF EPIPOLAR LINES FOR EACH REFERENCE POINT
	el1 <- c('l1x','l1y')
	el2 <- c('l2x','l2y')
	e_lines <- matrix(NA, nrow=nrow(curve_ref), ncol=ncol(curve_ref)+ncol(curve_1), 
		dimnames=list(NULL, c(el1, el2)))

	# GET EPIPOLAR LINES FROM REFERENCE VIEW IN OTHER VIEW
	for(i in 2:(nrow(curve_ref)-1)){

		# GET EPIPOLAR LINES
		e_line <- dltEpipolarLine(p=curve_ref[i, ], cal.coeff1=cal_coeff)

		# ADD TO MATRIX
		e_lines[i, ] <- c(e_line$l1, e_line$l2)
	}

	# PROPORTIONAL SIZE OF WINDOW FOR FINDING CLOSEST POINT TO EPIPOLAR
	epi.win.size.frac <- 0.1

	# SET INITIAL WINDOW SIZE
	epi_win_size <- c(-round(epi.win.size.frac*nrow(curve_1))/2, round(epi.win.size.frac*nrow(curve_1))/2)

	# SET MIN WINDOW SIZES
	epi.win.min.size <- c(-round((epi.win.size.frac/2)*nrow(curve_1))/2, round((epi.win.size.frac/2)*nrow(curve_1))/2)
	
	match_idx <- rep(NA, nrow(curve_ref))
	match_error <- rep(NA, nrow(curve_ref))
	rev_match_idx <- rep(NA, nrow(curve_1))
	
	# ASSUME START AND END POINTS ARE MATCHING
	match_idx[1] <- 1
	match_idx[length(match_idx)] <- nrow(curve_1)
	rev_match_idx[1] <- 1
	rev_match_idx[length(rev_match_idx)] <- nrow(curve_ref)
	
	# INITIAL OTHER START POSITION
	prev_oth_rel_pos <- 1

	# FOR EACH POINT ON THE REFERENCE CURVE, FIND CLOSEST POINT TO EPIPOLAR ON OTHER CURVE
	
	# INITIAL POINT MATCHING INDEX
	n <- 0
	i <- 2
	
	while(i < (nrow(curve_ref)-1) && n < nrow(curve_ref)*3){
		
		# SKIP POINTS WHERE EPIPOLAR ANGLE IS LESS THAN DIRECT MATCH THRESHOLD
		if(tan_angles[i] < min.direct.tangency*(pi/180)){
			i <- i + 1
			n <- n + 1
			next
		}
		
		# FIND RELATIVE POSITION IN OTHER CURVE
		oth_rel_pos <- round(max((i / nrow(curve_ref))*nrow(curve_1), 1))
		
		# IF PREDICTED CORRESPONDING INDEX (WINDOW CENTER) IS LESS THAN PREVIOUS, MAKE SAME AS PREVIOUS
		if(oth_rel_pos < prev_oth_rel_pos) oth_rel_pos <- prev_oth_rel_pos

		# SET WINDOW
		epi_win <- c(oth_rel_pos + floor(epi_win_size[1]), oth_rel_pos + floor(epi_win_size[2]))
		
		# MAKE SURE LIMITS ARE WITHIN MATRIX LIMITS
		if(epi_win[1] < 1) epi_win[1] <- 1
		if(epi_win[2] > nrow(curve_1)) epi_win[2] <- nrow(curve_1)

		# FIND DISTANCES TO EPIPOLAR LINE
		dpl <- distancePointToLine(p=curve_1[epi_win[1]:epi_win[2], ], e_lines[i, el1], e_lines[i, el2])

		# CHECK IF TWO MINIMA ARE PRESENT
		if(diff(range(order(dpl)[1:(min(15,length(dpl)))])) > 30){

			# FIND FIRST MIN
			mins <- which.min(dpl)

			# SET NEIGHBORING VALUES TO NA
			dpl_blocked <- dpl
			dpl_blocked[max(1, mins[1]-10):min(mins[1]+10, length(dpl))] <- NA

			# FIND SECOND MIN
			mins[2] <- which.min(dpl_blocked)

			# COMPARE DISTANCE TO EXPECTED POSITION OF CORRESPONDING POINT
			min_idx <- mins[which.min(abs(mins+epi_win[1]-oth_rel_pos-1))]
			
			#cat("Two mins possible\n")
			#plot(cbind(epi_win[1]:epi_win[2], dpl), type='l')
			#abline(v=epi_win[1]+min_idx-1, col='red')
			#abline(v=oth_rel_pos, col='green')

		}else{

			# FIND MINIMUM INDEX
			min_idx <- which.min(dpl)
		}
		
		# SET MATCH INDEX TEMPORARILY
		match_temp <- epi_win[1]+min_idx-1
		
		if(match_temp == 1){
			match_temp <- 2
			if(epi_win[1] > 1){i <- i + 1;n <- n + 1;next}
		}
		if(match_temp == nrow(curve_1)){
			match_temp <- nrow(curve_1) - 1
			if(epi_win[2] < nrow(curve_1)){i <- i + 1;n <- n + 1;next}
		}

		#cat(paste0('\t\t\t\ti:', i, '; win size:', floor(epi_win_size[1]), '-', floor(epi_win_size[2]), ' (', floor(sum(abs(epi_win_size))), '); win: ', epi_win[1], '-', epi_win[2], ' (', diff(range(epi_win)), '); match: ', match_temp, '; error: ', round(dpl[min_idx], 1), '\n'))

		# IF MATCH IS CLOSE TO WINDOW EDGE EXPAND WINDOW SIZE
		if(epi_win[2] - match_temp < 20 && epi_win[2] < nrow(curve_1)) epi_win_size[2] <- epi_win_size[2] + 10
		if(match_temp - epi_win[1] < 20 && epi_win[1] > 1) epi_win_size[1] <- epi_win_size[1] - 10

		# IF MATCH IS FAR FROM WINDOW EDGE CONTRACT WINDOW SIZE
		if(match_temp - epi_win[1] > 30 && epi_win[1] > 1) epi_win_size[1] <- epi_win_size[1] + 10
		if(epi_win[2] - match_temp > 30 && epi_win[2] < nrow(curve_1)) epi_win_size[2] <- epi_win_size[2] - 10
		
		# CHECK THAT WINDOW IS NOT LESS THAN MINIMUM SIZE
		epi_win_size[1] <- min(epi_win_size[1], epi.win.min.size[1])
		epi_win_size[2] <- max(epi_win_size[2], epi.win.min.size[2])

		# CHECK IF MINIMUM IS AT EDGE OF WINDOW
		if(match_temp == epi_win[1]){
			n <- n + 1
			next
		}

		# CHECK THAT MINIMUM DISTANCE IS NOT TOO LARGE (10 pixels)
		if(dpl[min_idx] > direct.min){
			i <- i + 1
			n <- n + 1
			next
		}

		# SAVE MATCH AND ERROR
		match_idx[i] <- match_temp
		rev_match_idx[match_temp] <- i
		match_error[i] <- dpl[min_idx]

		# SAVE OTHER CURVE RELATIVE INDEX FOR FUTURE ITERATIONS
		prev_oth_rel_pos <- oth_rel_pos

		i <- i + 1
		n <- n + 1
	}

	match_cts <- sum(!is.na(match_idx)) - 2

	#cat('\n')
	#cat(paste0('\t\t\t\tMatch error (1), mean: ', round(mean(match_error, na.rm=TRUE),2), ' min: ', round(min(match_error, na.rm=TRUE),2), ', max: ', round(max(match_error, na.rm=TRUE),2), '\n'))

	# FIND MATCH POINTS IN ORDER OF INCREASING UNCERTAINTY
	match_order <- order(tan_angles, decreasing=TRUE)

	# FILL IN GAPS
	n <- 0
	for(i in match_order){

		# SKIP FIRST AND LAST
		if(i == 1 || i == nrow(curve_ref)) next

		# SKIP IF ALREADY MATCHED
		if(!is.na(match_idx[i])) next

		# SKIP IF TANGENCY ANGLE IS LESS THAN MIN THRESHOLD
		if(tan_angles[i] < min.fill.tangency*(pi/180)) next
	
		#cat(i, '\n')
	
		# FIND PREVIOUS CLOSEST NON NA MATCH
		na_omit <- na.omit(match_idx[1:(i-1)])
		prev_non_na <- na_omit[length(na_omit)]

		# FIND NEXT CLOSEST NON NA MATCH
		next_non_na <- na.omit(match_idx[i:length(match_idx)])[1]

		# SET WINDOW
		epi_win <- c(prev_non_na, next_non_na)
	
		#cat('\t', epi_win[1], ',', epi_win[2], '\n')

		# FIND DISTANCES TO EPIPOLAR LINE
		dpl <- distancePointToLine(p=curve_1[epi_win[1]:epi_win[2], ], e_lines[i, el1], e_lines[i, el2])

		if(min(dpl, na.rm=TRUE) > 3) next

		# CHOOSE SUBSET FOR RECONSTRUCTION INCREASING THE THRESHOLD SEEMS TO HELP PREVENT DEVIATION FROM CURVE
		rec_sub_idx <- which(dpl < 6)
		#rec_sub_idx <- 1:length(dpl)
	
		# RESOLVED MATCHES IMMEDIATELY ON EACH SIDE
		if(epi_win[2]-epi_win[1] == 0){
			match_idx[i] <- epi_win[1]
			rev_match_idx[epi_win[1]] <- i
			match_error[i] <- dpl
			next
		}

		# GET 2D COORDINATES OF CLOSEST PREVIOUS AND NEXT MATCHED POINTS
		prev_next_2d <- cbind(curve_ref[c(rev_match_idx[prev_non_na], rev_match_idx[next_non_na]), ], 
			curve_1[c(prev_non_na, next_non_na), ])

		# RECONSTRUCT LINE REFERENCE POINTS
		line_ref_pts <- dltReconstruct(cal_coeff, prev_next_2d)$coor.3d
	
		if(length(rec_sub_idx) == 0) next
		if(length(rec_sub_idx) == 1){
			if(min(dpl, na.rm=TRUE) < 2){
				match_idx[i] <- rec_sub_idx+epi_win[1]-1
				rev_match_idx[rec_sub_idx+epi_win[1]-1] <- i
				match_error[i] <- min(dpl, na.rm=TRUE)
			}
			next
		}

		#cat('\trange(dpl): ', min(dpl), ', ', max(dpl), '\n')
		#cat('\tlength(rec_sub_idx): ', length(rec_sub_idx), '\n')

		# CREATE MATRIX OF POINTS TO RECONSTRUCT
		rec_sub_mat <- cbind(matrix(curve_ref[i, ], nrow=length(rec_sub_idx), ncol=2, byrow=TRUE), curve_1[rec_sub_idx+epi_win[1]-1, ])

		# RECONSTRUCT POINTS
		rec_sub_pts <- dltReconstruct(cal_coeff, rec_sub_mat)$coor.3d

		# FIND DISTANCE FROM RECONSTRUCTED POINTS TO REFERENCE LINE
		dpl_rec <- distancePointToLine(rec_sub_pts, line_ref_pts[1, ], line_ref_pts[2, ])

		# FIND COMPOSITE ERROR SCORE WITH TWO DISTANCE METRICS
		ttl_dpl <- epi.err.weight*((dpl[rec_sub_idx]-min(dpl[rec_sub_idx], na.rm=TRUE))/(diff(range(dpl[rec_sub_idx], na.rm=TRUE)))) + 
			rec.err.weight*((dpl_rec-min(dpl_rec, na.rm=TRUE))/(diff(range(dpl_rec, na.rm=TRUE))))

		# FIND MINIMUM INDEX
		#min_idx <- rec_sub_idx[which.min(dpl_rec)[1]] + epi_win[1] - 1
		min_idx <- rec_sub_idx[which.min(ttl_dpl)[1]] + epi_win[1] - 1
	
		match_idx[i] <- min_idx
		rev_match_idx[min_idx] <- i
		#match_error[i] <- dpl[rec_sub_idx[which.min(dpl_rec)[1]]]
		match_error[i] <- dpl[rec_sub_idx[which.min(ttl_dpl)[1]]]
	
		#cat('\t', which.min(dpl)+epi_win[1]-1, ',', min_idx, ',', match_error[i], '\n')
	
		n <- n + 1

		#if(n > 100) break
	}

	# NUMBER MATCHED IN FILLING GAPS
	match_cts[2] <- sum(!is.na(match_idx)) - match_cts[1]

	#cat(paste0('\t\t\t\tMatch error (2), mean: ', round(mean(match_error, na.rm=TRUE),2), ' min: ', round(min(match_error, na.rm=TRUE),2), ', max: ', round(max(match_error, na.rm=TRUE),2), '\n'))

	# ASSIGN CORRESPONDING POINTS
	for(i in 1:nrow(curve_ref)) if(!is.na(match_idx[i])) curve_1_hom[i, ] <- curve_1[match_idx[i], ]

	if(FALSE){

		# PLOT BOTH CURVES, CHOOSE LAYOUT BASED ON ASPECT RATIO
		if(diff(range(curve_ref[, 1])) / diff(range(curve_ref[1, ]))){layout(cbind(1,2))}else{layout(rbind(1,2))}

		par(mar=c(0,0,0,0))
		
		rainbow_grad <- rainbow(nrow(curve_ref), end=0.85)
		rainbow_alpha_grad <- rainbow(nrow(curve_ref), end=0.85, alpha=0.4)
		
		# PLOT REFERENCE CURVE
		plot(curve_ref, type='n', asp=1, yaxt='n', ylab='', xlab='')
		points(curve_ref, cex=0.1, col=rainbow_grad)
		points(curve_ref[c(1,nrow(curve_ref)),], cex=1, col=rainbow_grad[c(1, length(rainbow_grad))])

		# CREATE PLOT SPACE FOR OTHER CURVE
		plot(curve_1, type='n', asp=1, xaxt='n', yaxt='n', ylab='', xlab='')

		# PLOT
		if(FALSE){

			epi_x_range <- diff(range(curve_ref))*0.04

			for(i in seq(from=1, to=nrow(curve_ref), by=round(nrow(curve_ref)*0.01))){

				# FIND EQUATION FOR EPIPOLAR LINE
				m <- (e_lines[i, 'l2y'] - e_lines[i, 'l1y']) / (e_lines[i, 'l2x'] - e_lines[i, 'l1x'])
				b <- e_lines[i, 'l2y'] - m*e_lines[i, 'l2x']
			
				if(is.na(curve_1_hom[i, 1])) next

				#
				t <- seq(curve_1_hom[i, 1]-epi_x_range, curve_1_hom[i, 1]+epi_x_range, length=2)
				points(cbind(t, m*t + b), type='l', col=rainbow_alpha_grad[i])
			
			}
		}

		#
		error_cex <- (match_error - min(match_error, na.rm=TRUE)) / diff(range(match_error, na.rm=TRUE)) + 0.2

		# PLOT POINTS
		points(curve_1[c(1,nrow(curve_1)),], cex=1, col=rainbow_grad[c(1, length(rainbow_grad))])
		points(curve_1, cex=0.1, col=rgb(0,0,0,0.1))
		points(curve_1, cex=error_cex, pch=20, col=rainbow_grad[rev_match_idx])
	}
	
	# NUMBER UNMATCHED
	match_cts[3] <- sum(is.na(match_idx))
	
	names(match_cts) <- c('Initial match', 'Fill match', 'Unmatched')

	list('match'=curve_1_hom, 'match.cts'=match_cts)
}
