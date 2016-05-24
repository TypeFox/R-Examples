readLandmarksToList <- function(file, semilandmark.pattern='[0-9]+$', ...){

	# TO GROUP SEMILANDMARKS INCLUDE A PATTERN TO IDENTIFY ENUMERATION SCHEME
	#	DEFAULT IS A SERIES OF NUMBERS AT THE END OF THE LANDMARK NAME
	# TO NOT GROUP SEMILANDMARKS AND TREAT ALL AS SIMPLE LANDMARKS INPUT A BLANK PATTERN ('') OR FALSE

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

	# GET UNIQUE LANDMARK NAMES WITHOUT TERMINAL NUMBERS
	landmark_names <- rep(NA, 0)
	for(i in 1:length(landmark_file_list)) landmark_names <- c(landmark_names, unique(sub(pattern=semilandmark.pattern, replacement='', x=rownames(landmark_file_list[[i]]))))
	landmark_names <- sort(unique(landmark_names))

	# CREATE LANDMARK LIST WITH LANDMARK NAMES
	for(landmark_name in landmark_names) landmark_list[[landmark_name]] <- list()

	# SAVE LANDMARKS TO LANDMARK LIST
	for(i in 1:length(landmark_file_list)){
		for(j in 1:length(landmark_list)){
			
			# GET ROWS CORRESPONDING TO LANDMARK
			landmark_idx <- sub(pattern=semilandmark.pattern, replacement='', x=rownames(landmark_file_list[[i]])) == names(landmark_list)[j]

			# LANDMARK(S) NOT PRESENT IN CURRENT FILE
			if(sum(landmark_idx) == 0) next

			# ADD MATRIX TO LANDMARK LIST
			landmark_list[[names(landmark_list)[j]]][[i]] <- landmark_file_list[[i]][landmark_idx, ]
			
			# SKIP IF NOT CURVE POINTS
			if(is.null(nrow(landmark_list[[names(landmark_list)[j]]][[i]]))) next
			
			# GET ROWNAMES
			row_names <- rownames(landmark_list[[names(landmark_list)[j]]][[i]])

			# GET SEMILANDMARK NUMERIC SEQUENCE
			regexpr_r <- regexpr(pattern=semilandmark.pattern, text=row_names)

			# ISOLATE SEMILANDMARK NUMERIC SEQUENCE
			curve_pt_num <- as.numeric(substr(row_names, regexpr_r, regexpr_r+attr(regexpr_r, 'match.length')))
			
			# PAIR INDICES AND SEQUENCES FOR SORTING
			idx_seq <- matrix(1:length(row_names), nrow=length(row_names), ncol=1, dimnames=list(curve_pt_num, NULL))

			# SORT CURVE POINTS BY SEQUENCE
			landmark_list[[names(landmark_list)[j]]][[i]] <- landmark_list[[names(landmark_list)[j]]][[i]][idx_seq[as.character(sort(curve_pt_num)), ], ]
		}
	}
	
	landmark_list
}