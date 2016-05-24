copyTransformation <- function(m1, m2, mn, translate = TRUE, lar.cons = NULL, lar.compare = NULL){

	# IF ONLY ONE POINT PROVIDED ONLY APPLY TRANSLATION
	if(!is.matrix(m1)){
		if(is.matrix(mn)){
			return(mn + matrix(m2 - m1, nrow=nrow(mn), ncol=ncol(mn), byrow=TRUE))
		}else{
			return(mn + (m2 - m1))
		}
	}else{
		if(nrow(m1) != nrow(m2)) stop(paste0("'m1' ", nrow(m1), " must have the same number of rows as 'm2' ", nrow(m2), "."))
	}

	# CONVERT TO MATRIX IF VECTOR
	if(is.vector(mn)) mn <- matrix(mn, 1, 3)

	# IF AT LEAST TWO POINTS IN M1 ARE THE SAME, ONLY APPLY TRANSLATION
	if(sum(distPointToPoint(m1, m1[c(2:nrow(m1),1), ]) == 0) > 0)
		return(mn + matrix(m2[1, ] - m1[1, ], nrow=nrow(mn), ncol=ncol(mn), byrow=TRUE))

	# EMPTY RETURN MATRIX
	if(nrow(mn) == 0) return(mn)

	# TRANSFORMED POINT MATRIX
	mr <- matrix(NA, nrow=nrow(mn), ncol=ncol(mn))

	## CHECK WHETHER SIMPLE TRANSLATION COPIES TRANSFORMATION
	if(sum(apply(m2-m1, 2, 'sd')) < 1e-10){
		return(mn + matrix((m2-m1)[1, ], nrow=nrow(mn), ncol=3, byrow=TRUE))
	}
	

	# TRANSFORMATION GIVEN ONLY TWO POINTS (MINIMIZE ROTATION ABOUT LONG AXIS)
	two_ref_pts <- ifelse(nrow(m1) == 2, TRUE, FALSE)

	if(two_ref_pts){
		m1 <- rbind(m1, m1[1, ] + uvector(vorthogonal(m1[2, ] - m1[1, ])))
		m2 <- rbind(m2, m2[1, ] + uvector(vorthogonal(m2[2, ] - m2[1, ])))
	}

	# INCLUDE ORIGIN
	if(!translate) mn <- rbind(mn, c(0,0,0))

	# FIND INITIAL COORDINATE SYSTEM AXES
	rm_i <- setCoordinateAxes(m1)

	# FIND INITIAL POSITION VECTOR MATRIX
	mn_ref <- (mn - matrix(m1[1, ], nrow=nrow(mn), ncol=3, byrow=TRUE)) %*% t(rm_i)

	# FIND FINAL COORDINATE SYSTEM AXES
	ca <- setCoordinateAxes(m2)

	# FIND ROTATION MATRIX TO TRANSFORM TO FINAL COORDINATE SYSTEM
	rm <- tMatrixDC(ca)
	
	# TRANSFORM POINTS
	mr <- (mn_ref %*% t(rm)) + matrix(m2[1, ], nrow=nrow(mn_ref), ncol=3, byrow=TRUE)

	# CENTER BY PREVIOUS ORIGIN
	if(!translate){
		mr <- mr - matrix(mr[nrow(mr), ], nrow=nrow(mr), ncol=3, byrow=TRUE)
		mr <- mr[1:(nrow(mr)-1), ]
	}

	if(two_ref_pts){
		
		if(!is.null(lar.cons)){

			# TRANSFORM REFERENCE POINT
			pt_ref <- (((lar.cons$point - m1[1, ]) %*% t(rm_i)) %*% t(rm)) + m2[1, ]

			# FIND ROTATION ABOUT LONG AXIS THAT KEEPS CONSTRAINT
			# DEFINE CIRCLE USING POINT FOR RADIUS AND TWO JOINTS AS LONG-AXIS
			circle <- defineCircle(center=m2[1, ], nvector=m2[1, ]-m2[2, ], point_on_radius=pt_ref)

			# FIND INTERSECTION OF CIRCLE AND PLANE
			icp_t <- intersectCirclePlane(circle, P=lar.cons$point.i, N=lar.cons$vec)

			# ADD NEGATIVE ANGLES
			icp_t <- c(icp_t, -icp_t)

			# FIND THE DIRECTION IN WHICH TO ROTATE POINTS
			ref_rot <- matrix(NA, nrow=length(icp_t), ncol=3)
			for(i in 1:length(icp_t)) ref_rot[i, ] <- (pt_ref - m2[1, ]) %*% tMatrixEP(m2[1, ]-m2[2, ], a=icp_t[i]) + m2[1, ]
			dist_plane <- abs(distPointToPlane(ref_rot, n=lar.cons$vec, q=lar.cons$point.i))
			
			# KEEP PROPER ROTATION ANGLES
			icp_t <- icp_t[dist_plane < 1e-9]

			if(length(icp_t) > 1){

				# ROTATE POINTS OVER ANGLES
				ref_rot <- matrix(NA, nrow=length(icp_t), ncol=3)
				for(i in 1:length(icp_t)) ref_rot[i, ] <- (pt_ref - m2[1, ]) %*% tMatrixEP(m2[1, ]-m2[2, ], a=icp_t[i]) + m2[1, ]
				
				# FIND THE ANGLE THAT PUTS POINT CLOSEST TO THE PREVIOUS POINT
				icp_t <- icp_t[which.min(distPointToPoint(ref_rot, lar.compare))]
			}

			# CENTER POINTS ABOUT ONE OF POINTS ON AXIS OF ROTATION
			mr_centroid_m <- matrix(m2[1, ], nrow=nrow(mr), ncol=3, byrow=TRUE)
			mrc <- mr - mr_centroid_m

			# APPLY ROTATION
			mrr <- mrc %*% tMatrixEP(m2[1, ]-m2[2, ], a=icp_t)

			# UNDO CENTERING TRANSLATION
			mrr <- mrr + mr_centroid_m

			return(mrr)
		}
		
		# FIND A POINT NOT COINCIDENT WITH END POINTS
		ref_idx <- NA
		for(i in 1:nrow(mn)) if(distPointToLine(mn[i, ], m1[1, ], m1[2, ]) > 1e-13){ref_idx <- i;break}
			
		# IF ALL POINTS ARE COINCIDENT WITH END POINTS RETURN TRANSFORMED POINTS
		if(is.na(ref_idx)) return(mr)

		# FIND POINT ON LINE BETWEEN END POINTS, AT NORMAL TO REF POINT
		ref_norm <- pointNormalOnLine(mn[ref_idx, ], m1[1, ], m1[2, ])
		#print(distPointToLine(ref_norm, m1[1, ], m1[2, ]))
		
		# GET DISTANCE FROM LINE TO REF POINT
		ref_norm_dist <- distPointToPoint(mn[ref_idx, ], ref_norm)

		# TRANSFORM REF NORMAL USING SAME TRANSFORMATION AS TO MN
		ref_norm_t <- (((ref_norm - m1[1, ]) %*% t(rm_i)) %*% t(rm)) + m2[1, ]
		
		# VECTOR FROM REF NORM TO REF POINT
		ref_vector <- uvector(mn[ref_idx, ] - ref_norm)*ref_norm_dist

		# TRANSLATED REF POINT FROM REF NORM
		ref_toproj <- ref_proj <- ref_norm_t + ref_vector
		
		# FIND MINIMUM DISTANCE BETWEEN POINT AND TRANSFORMED POINT POTENTIAL LOCATIONS
		minDist <- distPointToPlane(p=ref_toproj, n=uvector(m2[2, ] - m2[1, ]), q=ref_norm_t)
		
		# PROJECT TRANSLATED REF POINT INTO PLANE
		if(minDist > 1e-10){

			ref_proj <- pointPlaneProj(q=ref_toproj, p=ref_norm_t, n=uvector(m2[2, ] - m2[1, ]))			
			#ref_proj <- ref_norm_t + uvector(ref_proj - ref_norm_t)*ref_norm_dist
		}

		# FIND ANGLE BETWEEN TRANSFORMED POINT AND PROJECTED POINT
		a <- avectors(mr[ref_idx, ] - ref_norm_t, ref_proj - ref_norm_t)
		
		# IF DIFFERENCE IN ANGLE IS 0, RETURN TRANSFORMED MATRIX AS IS
		if(a == 0) return(mr)

		# CENTER POINTS ABOUT ONE OF POINTS ON AXIS OF ROTATION
		mr_centroid_m <- matrix(m2[2, ], nrow=nrow(mr), ncol=3, byrow=TRUE)
		mrc <- mr - mr_centroid_m

		# APPLY ROTATIONS IN BOTH DIRECTIONS
		mrc_p <- mrc %*% tMatrixEP(m2[2, ]-m2[1, ], a=a)
		mrc_n <- mrc %*% tMatrixEP(m2[2, ]-m2[1, ], a=-a)

		# UNDO CENTERING TRANSLATION
		mrc_p <- mrc_p + mr_centroid_m
		mrc_n <- mrc_n + mr_centroid_m

		# IDENTIFY CORRECT DIRECTION OF ROTATION AND RETURN CORRESPONDING MATRIX
		ifelse(distPointToPoint(ref_proj, mrc_p[ref_idx, ]) < distPointToPoint(ref_proj, mrc_n[ref_idx, ]), return(mrc_p), return(mrc_n))
	}

	mr
}