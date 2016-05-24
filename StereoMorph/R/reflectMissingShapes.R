reflectMissingShapes <- function(shapes, file = NULL, print.progress = FALSE, ...){

	if(!is.list(shapes)){

		# FIND INPUT FILE PATHS		
		if(grepl('[.]txt$', shapes[1], ignore.case=TRUE)){
			shapes_in_fpaths <- shapes
		}else{
		
			if(is.null(file)) stop("If 'shapes' is a directory, file must also be a directory.")
		
			shapes_in_fpaths <- paste0(gsub('/$', '', shapes), '/', list.files(shapes))
		}
		
		if(is.null(file)){
			rms <- reflectMissingShapes(shapes=XML4R2list(shapes_in_fpaths[1]), file=NULL, print.progress=print.progress, ...)
			return(rms)
		}

		# SET OUTPUT FILE PATHS		
		if(grepl('[.]txt$', file[1], ignore.case=TRUE)){
			shapes_save_fpaths <- file
		}else{
			if(!file.exists(file)) dir.create(file)
			shapes_save_fpaths <- paste0(gsub('/$', '', file), '/', list.files(shapes))
		}

		# CALL DRAW SHAPES FOR EACH FILE
		for(i in 1:length(shapes_in_fpaths))
			rms <- reflectMissingShapes(shapes=XML4R2list(shapes_in_fpaths[i]), file=shapes_save_fpaths[i], print.progress=print.progress, ...)

		if(length(shapes_in_fpaths) == 1) return(rms)

		return(NULL)
	}

	if(length(file) > 1) stop("If 'shapes' input is shapes structure, file must be a vector of length 1 or NULL.")
	
	if('shapes' %in% names(shapes)){
		shapes_in <- shapes$shapes
	}else{
		shapes_in <- shapes
	}	

	# CREATE EMPTY CURVES LIST
	curves <- list()
	
	# CREATE MATRIX WITH ALL LANDMARKS AND CURVE POINTS
	lm.matrix <- shapes_in$landmarks

	# ADD CURVE POINTS
	if(!is.null(shapes_in$curves) && length(shapes_in$curves) > 0){
		for(i in 1:length(shapes_in$curves)){

			# ADD CURVE POINT NAMES
			rownames(shapes_in$curves[[i]]) <- paste0(names(shapes_in$curves)[i], '', formatC(1:nrow(shapes_in$curves[[i]]), width=6, format="d", flag="0"))
			
			# ADD TO MATRIX
			lm.matrix <- rbind(lm.matrix, shapes_in$curves[[i]])
		}
	}
	
	# REFLECT MISSING LANDMARKS AND CURVE POINTS
	reflect_missing <- reflectMissingLandmarks(lm.matrix, ...)
	
	if(print.progress){
		cat('\nreflectMissingShapes\n')
		if(!is.null(file)) cat('\tSave shapes as: ', file, '\n', sep='')
		
		cat('\tAlignment error:\n')
		cat('\t\t', paste(names(reflect_missing$align.error), ': ', format(reflect_missing$align.error), collapse='\n\t\t', sep=''), sep='')
		cat('\n')
	}

	lm.matrix <- reflect_missing$lm.matrix
	
	# FIND CURVE POINT ROW NAMES
	is_curve_pt_name <- grepl('[A-Za-z_][0-9]{6}$', rownames(lm.matrix))

	if(sum(is_curve_pt_name) > 0){
	
		# GET ALL CURVE POINT NAMES
		curve_point_names <- rownames(lm.matrix)[is_curve_pt_name]
		
		# EXTRACT CURVE NAME
		curve_name <- gsub('[0-9]{6}$', '', curve_point_names)

		# FIND UNIQUE CURVE NAMES
		curve_name_unique <- unique(curve_name)
		
		for(i in 1:length(curve_name_unique)){

			# GET CURVE POINTS
			curves[[curve_name_unique[i]]] <- lm.matrix[curve_point_names[curve_name == curve_name_unique[i]], ]
			
			# REMOVE ROWNAMES
			rownames(curves[[curve_name_unique[i]]]) <- NULL
			
			# REMOVE CURVE POINTS FROM LANDMARK MATRIX
			lm.matrix <- lm.matrix[!rownames(lm.matrix) %in% curve_point_names[curve_name == curve_name_unique[i]], ]
		}
	}

	# SHAPES LIST
	shapes <- list(
		'landmarks'=lm.matrix,
		'curves'=curves
		)

	if(!is.null(file)) list2XML4R(list('shapes'=shapes), file=file)
	
	class(shapes) <- 'shapes'

	return(shapes)
}