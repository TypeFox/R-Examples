landmarkListToMatrix <- function(lm.list){

	# GET UNIQUE LANDMARK NAMES
	landmark_names <- rep(NA, 0)
	num_sets <- rep(0, length(names(lm.list)))
	for(i in 1:length(names(lm.list))){

		if(is.null(lm.list[[names(lm.list)[i]]])) next
		
		# IF ELEMENT IS A MATRIX, MAKE FIRST AND ONLY ELEMENT
		if(is.matrix(lm.list[[names(lm.list)[i]]])) lm.list[[names(lm.list)[i]]] <- list(lm.list[[names(lm.list)[i]]])

		# GET NUMBER OF LANDMARK SETS
		num_sets[i] <- length(lm.list[[names(lm.list)[i]]])

		for(j in 1:length(lm.list[[names(lm.list)[i]]])){
			if(is.null(lm.list[[names(lm.list)[i]]][[j]])) next

			# IF SINGLE POINTS TAKE NAME OF SUB-LIST AND LENGTH AS NUMBER OF COLUMNS
			if(!is.matrix(lm.list[[names(lm.list)[i]]][[j]])){
				landmark_names <- c(landmark_names, names(lm.list)[i])
				num_col <- length(lm.list[[names(lm.list)[i]]][[j]])
				next
			}

			# GET NUMBER OF COLUMNS
			num_col <- ncol(lm.list[[names(lm.list)[i]]][[j]])

			# TAKE ROWNAMES (INCL SEMILANDMARK NUMBERING SCHEME)
			landmark_names <- c(landmark_names, rownames(lm.list[[names(lm.list)[i]]][[j]]))
		}
	}
	landmark_names <- sort(unique(landmark_names))

	# GET MAX NUMBER OF SETS
	num_sets <- max(num_sets)

	# GENERATE COLUMN NAMES
	col_names <- paste(c('x', 'y', 'z')[1:num_col], as.vector(matrix(1:num_sets, nrow=num_col, ncol=num_sets, byrow=TRUE)), sep="")

	# CREATE LANDMARK MATRIX WITH LANDMARK NAMES
	landmark_matrix <- matrix(NA, nrow=length(landmark_names), ncol=num_sets*num_col, dimnames=list(landmark_names, col_names))

	# SAVE LANDMARKS TO MATRIX
	for(i in 1:length(names(lm.list))){
		if(is.null(lm.list[[names(lm.list)[i]]])) next
		
		for(j in 1:length(lm.list[[names(lm.list)[i]]])){
			if(is.null(lm.list[[names(lm.list)[i]]][[j]])) next

			if(!is.matrix(lm.list[[names(lm.list)[i]]][[j]])){
				landmark_matrix[names(lm.list)[i], ((j-1)*num_col+1):((j-1)*num_col+num_col)] <- lm.list[[names(lm.list)[i]]][[j]]
			}else{
				landmark_matrix[rownames(lm.list[[names(lm.list)[i]]][[j]]), ((j-1)*num_col+1):((j-1)*num_col+num_col)] <- lm.list[[names(lm.list)[i]]][[j]]
			}			
		}
	}

	landmark_matrix
}