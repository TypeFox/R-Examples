removeOutlierCorners <- function(int_corners, nx, ny){

	plot <- FALSE

	# REMOVE EXTREME OUTLIER ITERATIONS
	n_max <- 4
	for(n in 1:n_max){

		# USE PRINCIPAL COMPONENTS TO FIT TWO ORTHOGONAL LINES TO POINTS
		xyCov <- cov(int_corners)
		eigenVectors <- eigen(xyCov)$vectors

		# FIND CENTER
		corner_center <- colMeans(int_corners)

		# GET POINTS TO DEFINE LINES
		t <- seq(min(int_corners[, 1])-mean(int_corners[, 1]), max(int_corners[, 1])-mean(int_corners[, 1]), len=2)
		pc1 <- cbind(t + mean(int_corners[, 1]), (eigenVectors[2,1]/eigenVectors[1,1])*t + mean(int_corners[, 2]))
		pc2 <- cbind(t + mean(int_corners[, 1]), (eigenVectors[2,2]/eigenVectors[1,2])*t + mean(int_corners[, 2]))

		# FIND THE DISTANCE FROM EACH INTERNAL CORNER TO PC AXES
		dist_pc1 <- distancePointToLine(int_corners, l1=c(pc1[1,1], pc1[1,2]), l2=c(pc1[2,1], pc1[2,2]))
		dist_pc2 <- distancePointToLine(int_corners, l1=c(pc2[1,1], pc2[1,2]), l2=c(pc2[2,1], pc2[2,2]))
	
		# FIND INTERNAL CORNERS ALONG PC1
		int_corners_pc1 <- int_corners[order(dist_pc1)[1:max(nx,ny)], ]

		# FIND CLOSEST POINT AMONG SAMPLE
		closest_point_pc1 <- rep(NA, nrow(int_corners_pc1))
		for(i in 1:nrow(int_corners_pc1))
			closest_point_pc1[i] <- min(distancePointToPoint(int_corners_pc1[i, ], int_corners_pc1[-i, ]))

		# FIND INTERNAL CORNERS ALONG PC1
		int_corners_pc2 <- int_corners[order(dist_pc2)[1:min(nx,ny)], ]

		# FIND CLOSEST POINT AMONG SAMPLE
		closest_point_pc2 <- rep(NA, nrow(int_corners_pc2))
		for(i in 1:nrow(int_corners_pc2))
			closest_point_pc2[i] <- min(distancePointToPoint(int_corners_pc2[i, ], int_corners_pc2[-i, ]))
	
		# ESTIMATE SQUARE SIZES IN EACH DIRECTION
		square_dim <- c(median(closest_point_pc1), median(closest_point_pc2))
		b1_scale <- square_dim[2]*0.5*min(nx,ny)
		b2_scale <- square_dim[1]*0.5*max(nx,ny)

		# SET BOUNDING LINES
		b11 <- rbind(pc1[1, ] + b1_scale*uvector_SM(pc2[1,]-corner_center), pc1[2, ] + b1_scale*uvector_SM(pc2[1,]-corner_center))
		b12 <- rbind(pc1[1, ] + b1_scale*-uvector_SM(pc2[1,]-corner_center), pc1[2, ] + b1_scale*-uvector_SM(pc2[1,]-corner_center))
		b21 <- rbind(pc2[1, ] + b2_scale*uvector_SM(pc1[1,]-corner_center), pc2[2, ] + b2_scale*uvector_SM(pc1[1,]-corner_center))
		b22 <- rbind(pc2[1, ] + b2_scale*-uvector_SM(pc1[1,]-corner_center), pc2[2, ] + b2_scale*-uvector_SM(pc1[1,]-corner_center))

		# FIND DISTANCE FROM PC1 RELATIVE TO CORRESPONDING BOUND
		dist_pc1_rel <- dist_pc1 / b1_scale
		dist_pc2_rel <- dist_pc2 / b2_scale

		# USE MAX RELATIVE DISTANCE FROM PC1 OR PC2 AS DISTANCE SCORE
		dist_score <- apply(cbind(dist_pc1_rel, dist_pc2_rel), 1, 'max')
	
		# 
		outliers <- which(dist_score / median(dist_score) > 3)
		if(n < n_max && length(outliers) > 0 && nrow(int_corners) - length(outliers) >= nx*ny){
			int_corners <- int_corners[-outliers, ]
			#print(outliers)
			#print(median(dist_score))
			#print(dist_score / median(dist_score))
		}else{
			break
		}
	}

	# PLOT PC REGRESSION LINES
	if(plot){
		#plot(int_corners, asp=1, col=gray((dist_score - min(dist_score)) / (diff(range(dist_score)))))
		plot(int_corners, asp=1)
		#plot(c(pc1[,1], pc2[,1]), c(pc1[,2], pc2[,2]), type='n', asp=1)
		points(pc1, type='l', col='red')
		points(x=corner_center[1], y=corner_center[2], col='brown')
		points(pc2, type='l', col='blue')
		points(b11, type='l', col='pink')
		points(b12, type='l', col='pink')
		points(b21, type='l', col='lightblue')
		points(b22, type='l', col='lightblue')
	}

	# ONCE GOOD ESTIMATE OF SQUARE DIMENSIONS, CHECK THAT EACH POINT HAS EXPECTED NUMBER OF CORNERS IN VICINITY
	if(nrow(int_corners) > nx*ny){

		remove <- rep(FALSE, nrow(int_corners))

		# START WITH FURTHEST POINTS
		n <- 1
		for(i in order(dist_score)[length(dist_score):1]){

			#if(nrow(int_corners)-sum(remove) == nx*ny) break

			dpp <- distancePointToPoint(int_corners[i, ], int_corners[-i, ])
			
			if(sum(dpp < max(square_dim)*2) < 2){remove[i] <- TRUE;next}
			if(sum(dpp < max(square_dim)*3) < 5){remove[i] <- TRUE;next}
			if(nx*ny > 10 && sum(dpp < max(square_dim)*4) < 10){remove[i] <- TRUE;next}

			n <- n + 1
			if(n > 20) break
		}

		if(plot) points(int_corners[remove, ], col='red', cex=1.1)
		
		int_corners <- int_corners[!remove, ]
		dist_score <- dist_score[!remove]
	}

	# GET RANK THRESHOLD JUST BELOW MAXIMUM DISTANCE SET
	dist_score_rank <- rank(dist_score)
	rank_thresh <- sort(dist_score_rank)[nx*ny]
	
	# ONLY RETAIN ROWS AT OR BELOW THRESHOLD
	corners_trim <- int_corners[dist_score_rank <= rank_thresh, ]

	if(plot) points(corners_trim, col='orange', cex=0.4)
	
	corners_trim
}