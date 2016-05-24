readLandmarksToArray <- function(file, na.omit=FALSE, ...){

	landmark_file_list <- list()
	file_list_is_null <- rep(NA, 0)
	landmark_list <- list()

	if(is.vector(file)){

		for(i in 1:length(file)){
			
			# SET DEFAULT NON-NULL LIST MATRIX
			file_list_is_null <- c(file_list_is_null, FALSE)

			# ERROR FOR NON-EXISTANT FILES
			if(!file.exists(file[i])) stop(paste0("File '", file[i], "' not found."))
	
			# ERROR FOR EMPTY FILES
			if(!file.info(file[i])$size) stop(paste0("File '", file[i], "' is empty."))
	
			# READ MATRIX
			read_table <- as.matrix(read.table(file[i], ...))
			
			# GET MATRIX DIMENSIONS
			mdim <- dim(read_table)

			# COPY LANDMARKS TO LIST
			landmark_file_list[[i]] <- read_table
		}
	}else{

		# TEMPORARY LIST FOR ARRAY FOR EACH ROW IN FILE PATH MATRIX
		landmark_array_list <- list()
		
		# SEND EACH ROW OF FILE MATRIX RECURSIVELY AND THEN ASSEMBLE INTO 4D ARRAY
		for(i in 1:nrow(file)){
			landmark_array_list[[i]] <- readLandmarksToArray(file[i, ], ...)
		}
		
		# ROWNAMES VECTOR
		row_names <- NULL

		# IS ARRAY EMPTY VECTOR
		is_array_empty <- rep(FALSE, length(landmark_array_list))
		
		for(i in 1:length(landmark_array_list)){

			# SKIP EMPTY ARRAYS
			if(0 %in% dim(landmark_array_list[[i]])){is_array_empty[i] <- TRUE; next}
			
			# GET COLUMN NAMES
			col_names <- dimnames(landmark_array_list[[i]])[[2]]

			# GET MATRIX DIMENSIONS
			mdim <- dim(landmark_array_list[[i]])[1:2]

			# IF ARRAY HAS ROWNAMES ADD ROWNAMES TO VECTOR
			if(!is.null(dimnames(landmark_array_list[[i]])[[1]])) row_names <- c(row_names, dimnames(landmark_array_list[[i]])[[1]])
		}

		if(is.null(row_names)){
			# CREATE EMPTY ARRAY
			landmark_array <- array(NA, dim=c(mdim[1], mdim[2], sum(!is_array_empty), ncol(file)), dimnames=list(NULL, col_names, NULL, NULL))
			
			j <- 1
			for(i in 1:length(landmark_array_list)){

				# SKIP EMPTY ARRAYS
				if(0 %in% dim(landmark_array_list[[i]])) next

				# FILL ARRAY
				landmark_array[, , j, ] <- landmark_array_list[[i]]
				
				j <- j + 1
			}
		}else{

			# CREATE EMPTY ARRAY
			landmark_array <- array(NA, dim=c(length(unique(row_names)), mdim[2], sum(!is_array_empty), ncol(file)), dimnames=list(sort(unique(row_names)), col_names, NULL, NULL))

			j <- 1
			for(i in 1:length(landmark_array_list)){

				# SKIP EMPTY ARRAYS
				if(0 %in% dim(landmark_array_list[[i]])) next

				# FILL ARRAY
				landmark_array[dimnames(landmark_array_list[[i]])[[1]], , j, ] <- landmark_array_list[[i]]
				
				j <- j + 1
			}
		}

		return(landmark_array)
	}

	# GET MATRIX COLUMN NUMBER (FIRST NON-NULL MATRIX IN LIST)
	for(i in 1:length(file)) if(!file_list_is_null[i]){mdim <- dim(landmark_file_list[[i]]);break}

	# GET UNIQUE LANDMARK NAMES
	landmark_names <- rep(NA, 0)
	for(i in 1:length(file)) if(!file_list_is_null[i]) landmark_names <- c(landmark_names, unique(rownames(landmark_file_list[[i]])))
	landmark_names <- sort(unique(landmark_names))

	if(!length(landmark_names)){

		# CREATE LANDMARK ARRAY ASSUMING EQUAL NUMBERS OF ROWS FOR EACH FILE
		landmark_array <- array(NA, dim=c(mdim[1], mdim[2], length(file)), dimnames=list(NULL, c('x', 'y', 'z')[1:mdim[2]], 1:length(file)))

		for(i in 1:length(file)){
			# SKIP NULL
			if(file_list_is_null[i]) next

			# SAVE LANDMARKS TO ARRAY
			landmark_array[, , i] <- landmark_file_list[[i]]
		}
	}else{
		# CREATE LANDMARK ARRAY WITH LANDMARK NAMES
		landmark_array <- array(NA, dim=c(length(landmark_names), mdim[2], length(file_list_is_null)), dimnames=list(landmark_names, c('x', 'y', 'z')[1:mdim[2]], 1:length(file)))
	
		for(i in 1:length(file)){
			# SKIP NULL
			if(file_list_is_null[i]) next

			# SAVE LANDMARKS TO ARRAY
			landmark_array[rownames(landmark_file_list[[i]]), , i] <- landmark_file_list[[i]]
		}
	}

	if(na.omit){
		if(dim(landmark_array)[3] == 1){
			non_na <- is.na(landmark_array[, 1, 1]) == 0
		}else{
			non_na <- rowSums(is.na(landmark_array[, 1, ])) == 0
		}

		# ERROR MESSAGE IF ALL ROWS HAVE AT LEAST ONE NA VALUE
		#if(sum(non_na) == 0) stop("With na.omit = TRUE return array is zero length, every row has at least one NA value.")

		# TRIM ARRAY TO ONLY ROWS WITHOUT ANY NA VALUES
		landmark_array <- landmark_array[non_na, , ]
	}

	landmark_array
}