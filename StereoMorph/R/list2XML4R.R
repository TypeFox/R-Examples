list2XML4R <- function(list, file="", ind=0){

	## CONVERTS LIST INTO XML4R STRUCTURE

	# XML4R STRING
	str <- ''

	if(length(list) == 0) return()

	for(i in 1:length(list)){

		name <- names(list)[i]

		if(is.null(list[[name]])) next

		if(typeof(list[[name]]) == 'list'){

			# OPENING TAG
			str <- c(str, paste0(paste(rep('\t', ind), collapse=''), '<', name, ' type=list>\n'))

			# SEND LIST RECURSIVELY
			str <- c(str, list2XML4R(list[[name]], ind=ind+1))
		}

		if(class(list[[name]]) == 'logical' || class(list[[name]]) == 'character' || class(list[[name]]) == 'integer' || class(list[[name]]) == 'vector' || class(list[[name]]) == 'numeric'){
			
			names <- ifelse(is.null(names(list[[name]])), FALSE, TRUE)
			
			type <- 'vector'
			if(class(list[[name]]) == 'logical') type <- 'logical'

			new_line <- '\n'
			if(length(list[[name]]) == 1) new_line <- ''

			# OPENING TAG
			str <- c(str, paste0(paste(rep('\t', ind), collapse=''), '<', name, ' type=', type, ' names=', names,' length=', length(list[[name]]),' as.numeric=', is.numeric(c(list[[name]])),' >', new_line))

			# VECTOR VALUES
			if(length(list[[name]]) == 1 && !names){
				str <- c(str, list[[name]])
			}else{
				if(names) str <- c(str, '\n', paste(rep('\t', ind+1), collapse=''), paste(names(list[[name]]), collapse='\t'), '\n')
				str <- c(str, paste(paste(rep('\t', ind), collapse=''), '\t', list[[name]], collapse='\n', sep=''), '\n')
			}

			# CLOSING TAG
			if(length(list[[name]]) == 1 && !names){
				str <- c(str, paste0('</', name, '>\n'))
				next
			}
		}

		if(class(list[[name]]) == 'matrix'){
		
			row_names <- ifelse(is.null(rownames(list[[name]])), FALSE, TRUE)
			col_names <- ifelse(is.null(colnames(list[[name]])), FALSE, TRUE)

			# OPENING TAG
			str <- c(str, paste0(paste(rep('\t', ind), collapse=''), '<', name, ' type=matrix rownames=', row_names,' colnames=', col_names,' nrow=', nrow(list[[name]]),' ncol=', ncol(list[[name]]),' as.numeric=', is.numeric(c(list[[name]])),' >\n'))

			if(col_names){
				str <- c(str, '\t', paste(rep('\t', ind), collapse=''))
				if(row_names) str <- c(str, '\t')
				str <- c(str, paste(colnames(list[[name]]), collapse='\t'), '\n')
			}

			# MATRIX VALUES
			for(r in 1:nrow(list[[name]])){
				str <- c(str, '\t', paste(rep('\t', ind), collapse=''))
				if(row_names) str <- c(str, rownames(list[[name]])[r], '\t')
				str <- c(str, paste(list[[name]][r, ], collapse='\t'), '\n')
			}
		}

		if(class(list[[name]]) == 'array'){
			
			# ARRAY DIMENSIONS
			dims <- dim(list[[name]])
			
			# CONVERT TO VECTOR
			v <- c(list[[name]])

			# OPENING TAG
			str <- c(str, paste0(paste(rep('\t', ind), collapse=''), '<', name, ' type=array dim=', paste(dims, collapse=','), ' as.numeric=', is.numeric(c(list[[name]])),' >\n'))
			
			if(length(dims) > 4) stop("listToXML4R() currently only works with arrays of 3 or 4 dimensions.")

			# ARRAY VALUES
			if(length(dims) == 3){
				for(j in 1:dims[3]){
					for(r in 1:dims[1]) str <- c(str, '\t', paste(rep('\t', ind), collapse=''), paste(list[[name]][r, , j], collapse='\t'), '\n')
					if(j < dims[3]) str <- c(str, '\n')
				}
			}

			# ARRAY VALUES
			if(length(dims) == 4){
				for(k in 1:dims[4]){
					for(j in 1:dims[3]){
						for(r in 1:dims[1]) str <- c(str, '\t', paste(rep('\t', ind), collapse=''), paste(list[[name]][r, , j, k], collapse='\t'), '\n')
						if(k < dims[4] || j < dims[3]) str <- c(str, '\n')
					}
				}
			}
		}

		# CLOSING TAG
		str <- c(str, paste0(paste(rep('\t', ind), collapse=''), '</', name, '>\n'))
	}

	str <- paste(str, collapse='')

	if(file != "") write(str, file)

	str
}
