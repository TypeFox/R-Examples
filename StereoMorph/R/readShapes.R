readShapes <- function(file, fields=NULL){
	
	# REPLACE EMTPY FILE NAMES WITH DIRECTORY
	file <- gsub('/.txt', '/', file)
	
	fields_internal <- fields
	if(!is.null(fields_internal)){

		# ADD FIELDS FOR INTERNAL OPERATIONS
		fields_internal <- unique(c(fields, 'image.id'))
	}
	
	# IF DIRECTORY
	if(!grepl('[.]txt$', file[1])){

		# CHECK THAT DIRECTORY EXISTS
		if(!file.exists(file)) stop(paste0("'", file, "' not found."))

		# LIST FILES, IF DIRECTORY
		file <- paste0(gsub('/+$', '', file), '/', list.files(file))
	}
	
	# GET FILENAMES
	str_split <- strsplit(file, '/')
	last_name <- rep(NA, length(file))
	last2_name <- rep(NA, length(file))
	for(i in 1:length(str_split)){
		last_name[i] <- tail(str_split[[i]], 1)
		if(length(str_split[[i]]) > 1) last2_name[i] <- str_split[[i]][length(str_split[[i]])-1]
	}
	
	if(length(unique(last_name)) == length(file)){
		filenames <- gsub('[.]txt$', '', last_name)
	}else if(length(unique(last2_name)) == length(file)){
		filenames <- last2_name
	}else{
		filenames <- 1:length(file)
	}
	
	if(length(file) == 1){

		# READ SINGLE FILE
		rlist <- readXML4R(file)$shapes

	}else{
	
		dims <- list()
		read_shapes <- list()
		
		# FIRST GET DIMENSIONS OF ALL MATRICES ACROSS FILES
		for(i in 1:length(file)){

			# READ FILE INTO LIST
			read_xml4r <- readXML4R(file[i])

			# SKIP IF FILE IS EMPTY (CONTAINS NO SUB-SHAPE ELEMENTS)
			if(is.null(read_xml4r$shapes)) next
			
			# GET SHAPES
			read_shapes[[i]] <- read_xml4r$shapes
			
			for(name in names(read_shapes[[i]])){
				
				# CREATE ENTRY IF ONE DOES NOT ALREADY EXIST
				if(!name %in% names(dims)) dims[[name]] <- list()

				# ONLY GET DIMENSIONS FROM MATRICES
				if(class(read_shapes[[i]][[name]]) %in% c('character', 'numeric')){
					dims[[name]] <- list('type' = 'vector', 'length' = length(file))
				}

				if(class(read_shapes[[i]][[name]]) %in% c('matrix')){

					# CREATE EMPTY FIELDS
					if(!'nrow' %in% names(dims[[name]])){
						dims[[name]] <- list(
							'type' = 'matrix', 
							'nrow' = rep(NA, length(file)), 'ncol' = rep(NA, length(file)),
							'rownames' = c(), 'colnames' = c()
						)
					}
				
					# FILL FIELDS
					dims[[name]][['nrow']][i] <- nrow(read_shapes[[i]][[name]])
					dims[[name]][['ncol']][i] <- ncol(read_shapes[[i]][[name]])
					if(!is.null(rownames(read_shapes[[i]][[name]]))) dims[[name]][['rownames']] <- unique(c(dims[[name]][['rownames']], rownames(read_shapes[[i]][[name]])))
					if(!is.null(colnames(read_shapes[[i]][[name]]))) dims[[name]][['colnames']] <- unique(c(dims[[name]][['colnames']], colnames(read_shapes[[i]][[name]])))
				}

				if(class(read_shapes[[i]][[name]]) %in% c('list')){
					dims[[name]] <- list('type' = 'list', 'length' = length(file))
				}
			}
		}
		
		# CREATE RETURN LIST
		rlist <- list()
		for(name in names(dims)){

			if(!is.null(fields_internal)) if(!name %in% fields_internal) next

			if(dims[[name]][['type']] == 'vector'){
				rlist[[name]] <- rep(NA, length(file))
				names(rlist[[name]]) <- filenames
			}
				
			if(dims[[name]][['type']] == 'matrix'){

				# GET NUMBER OF ROWS
				nrow <- max(dims[[name]][['nrow']], na.rm=TRUE)
				if(!is.null(dims[[name]][['rownames']])) nrow <- length(dims[[name]][['rownames']])

				# GET NUMBER OF COLUMNS
				ncol <- max(dims[[name]][['ncol']], na.rm=TRUE)
				if(!is.null(dims[[name]][['colnames']])) ncol <- length(dims[[name]][['colnames']])

				# CREATE ARRAY
				rlist[[name]] <- array(NA, dim=c(nrow, ncol, length(file)), dimnames=list(dims[[name]][['rownames']], dims[[name]][['colnames']], filenames))
			}

			if(dims[[name]][['type']] == 'list'){

				rlist[[name]] <- list()
			}
		}

		# FILL RETURN LIST
		for(i in 1:length(read_shapes)){
		
			for(name in names(read_shapes[[i]])){
			
				if(!is.null(fields_internal)) if(!name %in% fields_internal) next

				if(dims[[name]][['type']] == 'vector') rlist[[name]][i] <- read_shapes[[i]][[name]]
				
				if(dims[[name]][['type']] == 'matrix'){
					
					row_idx <- 1:nrow(read_shapes[[i]][[name]])
					if(!is.null(rownames(read_shapes[[i]][[name]]))) row_idx <- rownames(read_shapes[[i]][[name]])

					col_idx <- 1:ncol(read_shapes[[i]][[name]])
					if(!is.null(colnames(read_shapes[[i]][[name]]))) col_idx <- colnames(read_shapes[[i]][[name]])
					
					# COPY IN VALUES
					rlist[[name]][row_idx, col_idx, i] <- read_shapes[[i]][[name]][row_idx, col_idx]
				}

				if(dims[[name]][['type']] == 'list'){
					if(length(read_shapes[[i]][[name]]) == 0){
						rlist[[name]][[filenames[i]]] <- NULL
					}else{
						rlist[[name]][[filenames[i]]] <- read_shapes[[i]][[name]]
					}
				}
			}
		}
	}

	# REMOVE EMPTY LISTS
	#for(name in names(rlist)) if(is.list(rlist[[name]]) && length(rlist[[name]])) rlist[[name]] <- NULL

	# REMOVE OBJECTS NOT IN FIELDS
	if(!is.null(fields)) for(name in names(rlist)) if(!name %in% fields) rlist[[name]] <- NULL
	for(name in names(rlist)) if(length(rlist[[name]]) == 0) rlist[[name]] <- NULL

	class(rlist) <- 'shapes'
	rlist
}

print.shapes <- function(x, ...){

	r <- ''
	r <- c(r, '\nShapes\n')

	vector_limit <- 3
	vector_names <- c('image.id', 'scaling', 'scaling.units', 'ruler.pixel', 'ruler.interval', 
		'checkerboard.nx', 'checkerboard.ny', 'square.pixel', 'square.size')

	for(vector_name in vector_names){

		if(is.null(x[[vector_name]])) next
			
		vector_min <- min(vector_limit, length(x[[vector_name]]))
		
		r <- c(r, '\t', vector_name, ': ', paste(format(x[[vector_name]][1:vector_min]), collapse=', '))
		if(length(x[[vector_name]]) > vector_min) r <- c(r, ', ...')
		r <- c(r, '\n')
	}

	matrix_limit <- 3
	matrix_names <- c('landmarks', 'landmarks.pixel', 'landmarks.scaled', 'ruler.points', 'checker.pixel')

	for(matrix_name in matrix_names){

		if(is.null(x[[matrix_name]])) next

		r <- c(r, '\t', matrix_name, '\n')

		r <- c(r, '\t\tDimensions: ', paste(dim(x[[matrix_name]]), collapse=" x "), '\n')

		dimnames_1 <- dimnames(x[[matrix_name]])[[1]]
		if(!is.null(dimnames_1)){
			matrix_min <- min(matrix_limit, length(dimnames_1))
			r <- c(r, '\t\tRownames: ', paste(dimnames_1[1:matrix_min], collapse=", "))
			if(length(dimnames_1) > matrix_min) r <- c(r, ', ...')
			r <- c(r, '\n')
		}

		if(length(dim(x[[matrix_name]])) == 3){
			dimnames_3 <- dimnames(x[[matrix_name]])[[3]]

			if(!is.null(dimnames_3)){
				matrix_min <- min(matrix_limit, length(dimnames_3))
				r <- c(r, '\t\tMatrix names: ', paste(dimnames_3[1:matrix_min], collapse=", "))
				if(length(dimnames_3) > matrix_min) r <- c(r, ', ...')
			}
			r <- c(r, '\n')
		}
	}

	list_limit <- 3
	list_names <- c('curves.control', 'curves.pixel', 'curves.scaled', 'curves')
	for(list_name in list_names){

		if(is.null(x[[list_name]])) next

		r <- c(r, '\t', list_name, '\n')

		if(!is.null(names(x[[list_name]]))){
			
			is_any_list <- FALSE
			for(name in names(x[[list_name]])) if(is.list(x[[list_name]][[name]])) is_any_list <- TRUE
		
			if(is_any_list){

				j <- 1
				for(name in names(x[[list_name]])){
					r <- c(r, '\t\t', name, '\n')

					list_min <- min(list_limit, length(x[[list_name]][[name]]))
					r <- c(r, '\t\t\tList names: ', paste(names(x[[list_name]][[name]])[1:list_min], collapse=", "))
					if(length(x[[list_name]][[name]]) > list_min) r <- c(r, ', ...')
					r <- c(r, '\n')

					if(j == list_limit && length(x[[list_name]]) > list_limit){
						r <- c(r, '\t\t...\n')
						break
					}
					j <- j + 1
				}
			}else{
				list_min <- min(list_limit, length(x[[list_name]]))
				r <- c(r, '\t\tList names: ', paste(names(x[[list_name]])[1:list_min], collapse=", "))
				if(length(x[[list_name]]) > list_min) r <- c(r, ', ...')
				r <- c(r, '\n')
			}
		}
	}
	
	cat(r, sep='')
}