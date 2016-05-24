###############################################################
##checks whether data satisfy regular minimality/maximality####
###############################################################
check.regular <-
function(X, type = c("probability.different", "percent.same", "reg.minimal", "reg.maximal")){
	# X: the data
	# type: whether regular minimality or regular maximality is to be checked; "probability.different" and 
	#       "percent.same" are for datasets in the probability-different and percent-same formats, and imply
	#       regular minimality and regular maximality checks, respectively; "reg.minimal" and "reg.maximal"
	#       can be specified to force checking regular minimality and regular maximality, respectively,
	#       independent of the used dataset; the latter two values, in particular, are to be used for datasets
	#       that are properly in the general format; the default is "probability.different"

	type <- match.arg(type)
	if (type == "probability.different"){
		M <- check.data(X, format = "probability.different")
		check.reg.minimal <- TRUE
	}
	if (type == "percent.same"){
		M <- check.data(X, format = "percent.same")
		check.reg.minimal <- FALSE
	}
	if (type == "reg.minimal"){
		M <- check.data(X, format = "general")
		check.reg.minimal <- TRUE
	}
	if (type == "reg.maximal"){
		M <- check.data(X, format = "general")
		check.reg.minimal <- FALSE
	}
	if (is.null(M))
		invisible(NULL)
	else{
		PSE <- matrix(nrow = dim(M)[1], ncol = 2)  # matrix of points of subjective equality (PSEs) 
		M.transformed <- matrix(nrow = dim(M)[1], ncol = dim(M)[2])  # transformed matrix for which regular minimality/maximality 
	                                                                 # is satisfied in canonical form

		if (check.reg.minimal){  # checks regular minimality
			for (i in 1:dim(M)[1]){
				index.row <- which(M[i, ] == min(M[i, ]))
				if (length(index.row) > 1)
					return(message(paste("regular minimality is violated: more than one minimal entry found in row #", i, sep = "")))
				index.col <- which(M[, index.row] == min(M[, index.row]))
				if (length(index.col) > 1)
					return(message(paste("regular minimality is violated: more than one minimal entry found in column #", index.row, sep = "")))
				if (index.col != i)
					return(message(paste("regular minimality is violated: entry in row #", i," and column #", index.row, 
				                         " is minimal in row #", i, " but not in column #", index.row, sep = "")))
				PSE[i, ] <- cbind(i, index.row)
			}
		} else{  # checks regular maximality
			for (i in 1:dim(M)[1]){
				index.row <- which(M[i, ] == max(M[i, ]))
				if (length(index.row) > 1)
					return(message(paste("regular maximality is violated: more than one maximal entry found in row #", i, sep = "")))
				index.col <- which(M[, index.row] == max(M[, index.row]))
				if (length(index.col) > 1)
					return(message(paste("regular maximality is violated: more than one maximal entry found in column #", index.row, sep = "")))
				if (index.col != i)
					return(message(paste("regular maximality is violated: entry in row #", i," and column #", index.row,
				                         " is maximal in row #", i, " but not in column #", index.row, sep = "")))
			 PSE[i, ] <- cbind(i, index.row)
			}
		}
		# permutes columns of matrix M so that regular minimality/maximality is satisfied in the canonical form (canonical representation),
		# and assigns identical labels to all PSE pairs, here leaving intact the labeling of the rows but relabeling the columns with their
		# corresponding PSEs (canonical relabeling)
		for (i in 1:nrow(PSE)){  # canonical representation
			M.transformed[, PSE[i, 1]] <- M[, PSE[i, 2]]
		}
		dimnames(M.transformed) <- dimnames(M)  # canonical relabeling

		# data frame of PSEs with their common labels
		for (i in 1:nrow(PSE))  # replaces row and column numbers by their labels
			PSE[PSE == i] <- rownames(M)[i]  
		PSE <- data.frame(observation.area.1 = PSE[, 1], observation.area.2 = PSE[, 2], common.label = PSE[, 1], stringsAsFactors = FALSE)

		return(list(canonical.representation = M.transformed, canonical.transformation = PSE,
	                check = ifelse(check.reg.minimal, "regular minimality", "regular maximality"),
	                in.canonical.form = identical(PSE[, 1], PSE[, 2])))
	}
}
