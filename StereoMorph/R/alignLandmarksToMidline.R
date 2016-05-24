alignLandmarksToMidline <- function(lm.matrix, left = '(_l[_]?|_left[_]?)([0-9]*$)', 
	right = '(_r[_]?|_right[_]?)([0-9]*$)', left.remove = '\\2', right.remove = '\\2', 
	average = FALSE){
	# Modified from R function AMP() written by Annat Haber
	# The function uses midline landmarks and the average of bilateral landmarks to find the midline
	# Then aligns all landmarks to the midline

	# MAKE ID VECTOR
	id_side <- rep(NA, length(rownames(lm.matrix)))

	# ID EACH LANDMARK AS LEFT, RIGHT OR MIDLINE
	id_side[grepl(pattern=left, x=rownames(lm.matrix), ignore.case=TRUE)] <- 'L'
	id_side[grepl(pattern=right, x=rownames(lm.matrix), ignore.case=TRUE)] <- 'R'
	id_side[is.na(id_side)] <- 'M'

	# GET LIST OF LANDMARK NAMES WITHOUT SIDES
	landmark_names <- gsub(pattern=left, replacement=left.remove, x=rownames(lm.matrix), ignore.case=TRUE)
	landmark_names <- gsub(pattern=right, replacement=right.remove, x=landmark_names, ignore.case=TRUE)

	# GET UNIQUE LIST OF LANDMARK NAMES WITHOUT SIDES
	unique_landmark_names <- unique(landmark_names)

	# MIDLINE LANDMARK MATRIX
	midline_matrix <- matrix(NA, nrow=0, ncol=ncol(lm.matrix))

	for(landmark_name in unique_landmark_names){

		# IF LANDMARK IS LEFT OR RIGHT BUT NOT BILATERAL, SKIP LANDMARK
		if(sum(landmark_name == landmark_names) == 1 && !('M' %in% id_side[landmark_name == landmark_names])) next

		# GET MIDLINE AND BILATERAL LANDMARKS
		midline_landmarks <- lm.matrix[landmark_name == landmark_names, ]

		# SKIP IF ANY LANDMARKS ARE NA
		if(sum(is.na(midline_landmarks)) > 0) next

		# GET MEAN OF BILATERAL LANDMARKS
		if(is.matrix(midline_landmarks)) midline_landmarks <- colMeans(midline_landmarks)

		# CONVERT TO MATRIX AND ADD ROWNAME
		midline_landmarks <- matrix(midline_landmarks, nrow=1);rownames(midline_landmarks) <- landmark_name

		# ADD LANDMARK TO MIDLINE MATRIX
		midline_matrix <- rbind(midline_matrix, midline_landmarks)
	}

	# IF NO MIDLINE POINTS IN MATRIX, RETURN NULL
	if(nrow(midline_matrix) == 0) return(NULL)

	# FIND CENTROID OF MIDLINE LANDMARKS
	midline_matrix_centroid <- matrix(colMeans(midline_matrix), byrow=TRUE, nrow=nrow(lm.matrix), ncol=ncol(lm.matrix))

	# CENTER ALL LANDMARKS BY CENTROID OF MIDLINE POINTS
	lm.matrix <- lm.matrix - midline_matrix_centroid

	# FIND ROTATION MATRIX TO APPLY TO ENTIRE LANDMARK MATRIX
	SVD <- svd(midline_matrix - midline_matrix_centroid[1:nrow(midline_matrix), ])

	# ROTATE LANDMARKS TO ALIGN WITH MID-SAGITTAL PLANE
	lm.matrix <- lm.matrix %*% SVD$v
	
	if(average){
	
		# SET MIDLINE POINTS TO ZERO
		lm.matrix[id_side == 'M', 3] <- 0
		
		lm_matrix_avg <- lm.matrix

		# FIND CORRESPONDING BILATERAL LANDMARKS
		for(i in 1:nrow(lm.matrix)){

			match <- which(landmark_names[i] == landmark_names)

			if(length(match) > 1){

				# GET AVERAGE POSITION
				new_pos <- colMeans(cbind(lm.matrix[match[1:2], 1:2], abs(lm.matrix[match[1:2], 3])), na.rm=TRUE)
				lm_matrix_avg[i, ] <- c(new_pos[1:2], new_pos[3]*sign(lm.matrix[i, 3]))
			}
		}

		# AVERAGE BILATERAL LANDMARKS
		lm.matrix <- lm_matrix_avg
	}

	# FIND MIDLINE ERRORS
	midline_error <- rep(NA, 0)
	midline_names <- rep(NA, 0)

	for(landmark_name in rownames(midline_matrix)){

		if(sum(landmark_name == landmark_names) != 1) next

		# GET MIDLINE ERROR, SQUARE OF Z-COORDINATE DIFFERENCE FROM ZERO
		midline_error <- c(midline_error, sqrt((lm.matrix[landmark_name == landmark_names, ncol(lm.matrix)])^2))
		midline_names <- c(midline_names, landmark_name)
	}

	# ADD NAMES TO ERROR VECTOR
	names(midline_error) <- midline_names

	l <- list(lm.matrix=lm.matrix, midline.error=midline_error)
	class(l) <- 'alignLandmarksToMidline'
	return(l)
}

summary.alignLandmarksToMidline <- function(object, ...){
	r <- ''
	r <- c(r, '\nalignLandmarksToMidline Summary\n')

	r <- c(r, '\tAlignment Midline RMSE: ', sqrt(sum(object$midline.error)/length(object$midline.error)), '\n')
	for(i in 1:length(object$midline.error)) r <- c(r, '\t\t', names(object$midline.error)[i], ": ", object$midline.error[i], '\n')

	class(r) <- "summary.alignLandmarksToMidline"
	r
}

print.summary.alignLandmarksToMidline <- function(x, ...) cat(x, sep='')