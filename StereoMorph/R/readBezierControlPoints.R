readBezierControlPoints <- function(file, ndim = 2, ...){

	bezier_file_list <- list()
	bezier_list <- list()

	# COPY CONTROL POINTS TO LIST
	for(i in 1:length(file)){

		# ERROR FOR NON-EXISTANT FILES
		if(!file.exists(file[i])) stop(paste0("File '", file[i], "' not found."))
	
		# ERROR FOR EMPTY FILES
		if(!file.info(file[i])$size) return(list())
	
		# FIND MAX NUMBER OF COLUMNS IN FILE
		Y <- readLines(file[i], warn = FALSE, ...)
		col_ct <- unlist(lapply(Y, function(x) length(unlist(strsplit(x, "\t")))))
	
		# SET NUMBER OF ROWS
		n_rows <- sum(col_ct > 0)
	
		# MAKE EMPTY MATRIX
		m <- matrix(NA, nrow=n_rows, ncol=max(col_ct))
	
		j <- 1
		for(line in Y){
	
			# GET ROW VALUES
			row <- strsplit(line, "\t")[[1]]
	
			# SKIP EMPTY ROWS
			if(length(row) == 0) next
			
			# REPLACE EMPTY CELLS WITH NA
			row[row == ""] <- NA
	
			# ADD ROW TO MATRIX
			m[j, 1:length(row)] <- row
			
			j <- j + 1
		}
	
		# REMOVE ROWS WITH ONLY NA
		m <- matrix(m[rowSums(!is.na(m)) > 0, ], ncol=ncol(m))
		
		# SEPARATE ROW NAMES FROM VALUES
		bezier_file_list[[i]] <- matrix(m[, 2:ncol(m)], nrow=nrow(m), ncol=ncol(m)-1, dimnames=list(m[, 1], NULL))
	}

	# CREATE CURVE LIST WITH CURVE NAMES
	for(i in 1:length(bezier_file_list)) for(curve_name in rownames(bezier_file_list[[i]])) bezier_list[[curve_name]] <- list()

	# SAVE CURVES TO CURVE LIST
	for(i in 1:length(bezier_file_list)) for(j in 1:nrow(bezier_file_list[[i]])){

		# CONVERT SINGLE ROW OF CONTROL POINTS TO NDIM-DIMENSIONAL MATRIX
		curve_m <- matrix(as.numeric(bezier_file_list[[i]][j, !is.na(bezier_file_list[[i]][j, ])]), ncol=ndim, byrow=T)

		# ADD MATRIX TO CURVE LIST
		bezier_list[[rownames(bezier_file_list[[i]])[j]]][[i]] <- curve_m
	}

	bezier_list
}