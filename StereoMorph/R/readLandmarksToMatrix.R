readLandmarksToMatrix <- function(file, na.omit = FALSE, ...){

	landmark_file_list <- list()
	landmark_list <- list()

	for(i in 1:length(file)){

		# ERROR FOR NON-EXISTANT FILES
		if(!file.exists(file[i])) stop(paste0("File '", file[i], "' not found."))

		# ERROR FOR EMPTY FILES
		if(!file.info(file[i])$size) stop(paste0("File '", file[i], "' is empty."))

		# COPY LANDMARKS TO LIST
		landmark_file_list[[i]] <- as.matrix(read.table(file[i], ...))
	}

	# GET MATRIX COLUMN NUMBER (FIRST NON-NULL MATRIX IN LIST)
	for(i in 1:length(landmark_file_list)) if(!is.null(landmark_file_list[[i]])){num_col <- ncol(landmark_file_list[[i]]);break}

	# GET UNIQUE LANDMARK NAMES
	landmark_names <- rep(NA, 0)
	for(i in 1:length(landmark_file_list)) landmark_names <- c(landmark_names, unique(rownames(landmark_file_list[[i]])))
	landmark_names <- sort(unique(landmark_names))

	# GENERATE COLUMN NAMES
	col_names <- paste(c('x', 'y', 'z')[1:num_col], as.vector(matrix(1:length(file), nrow=num_col, ncol=length(file), byrow=TRUE)), sep="")

	# CREATE LANDMARK MATRIX WITH LANDMARK NAMES
	landmark_matrix <- matrix(NA, nrow=length(landmark_names), ncol=length(file)*num_col, dimnames=list(landmark_names, col_names))

	# SAVE LANDMARKS TO MATRIX
	for(i in 1:length(landmark_file_list)) landmark_matrix[rownames(landmark_file_list[[i]]), ((i-1)*num_col+1):((i-1)*num_col+num_col)] <- landmark_file_list[[i]]

	# OMIT NA IF SPECIFIED
	if(na.omit){
		non_na <- rep(TRUE, nrow(landmark_matrix))

		# FIND ROWS THAT HAVE AT LEAST ONE NA
		for(i in 1:nrow(landmark_matrix)) if(sum(is.na(landmark_matrix[i, ])) > 0) non_na[i] <- FALSE

		# TRIM ARRAY TO ONLY ROWS WITHOUT ANY NA VALUES
		landmark_matrix <- landmark_matrix[non_na, ]
	}

	landmark_matrix
}