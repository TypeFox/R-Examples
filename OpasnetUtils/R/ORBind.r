# ORBIND ########################### orbind combines two ovariables using clever rbind

orbind <- function(x, y) {
	if(class(x) == "ovariable") x <- x@output
	if(class(y) == "ovariable") y <- y@output
	
	temp1 <- addmissingcol(x, y) # See below
	temp2 <- addmissingcol(y, x)
	
	return(rbind(temp1, temp2))
	
#Should this be made S4 function for ovariables? Then it could be named simply rbind.
	
}

addmissingcol <- function(e1, e2) { #Adds all missing columns. Merges Iter if that is missing.
	cols <- setdiff(colnames(e2), colnames(e1)) # Take all columns that do not exist in e1.
	
	if("Iter" %in% cols) {
		e1 <- merge(unique(e2["Iter"]), e1) # Add Iter with all locations existing in e2.
		cols <- cols[cols != "Iter"] # Remove Iter from the list of columns to add.
	}
	col <- as.data.frame(array(NA, dim = c(1, length(cols))))
	colnames(col) <- cols
	if("Unit" %in% cols) {col[, "Unit"] <- "?"}
	
	return(cbind(e1, col))
}

##############
# Combine 
################
# Combine ovariables, similar to orbind but has a different operating principle and allows multiple ovariables
# at once. Source is preserved but it is no longer a full marginal. 
# ... - any number of ovariables separated by commas
# name - of the resulting ovariable
#############


combine <- function(..., name = character()) {
	variable_list <- list(...)
	new_source_col <- paste(name, "Source", sep = "")
	new_res_col <- paste(name, "Result", sep = "")
	marginals <- character()
	for (i in 1:length(variable_list)) {
		
		var <- variable_list[[i]]
		if (class(var) != "ovariable") stop(paste("Variable #", i, "not ovariable."))
		if (nrow(var@output) == 0) stop(paste(var@name, "not evaluated."))
		old_source_col <- paste(var@name, "Source", sep = "")
		
		# Get index column names (excluding own source)
		temp <- colnames(var@output)[var@marginal & colnames(var@output) != old_source_col]
		marginals <- c(marginals, as.character(temp))
		
		# Reformat source columns
		temp <- var@output[[old_source_col]]
		colnames(var@output)[colnames(var@output) == old_source_col] <- new_source_col
		var@output[[new_source_col]] <- paste(var@name, temp, sep = "") # XFormula, XData, X, etc.
		
		# Rename individual result columns to common name
		
		colnames(var@output)[colnames(var@output) == paste(var@name, "Result", sep = "")] <- new_res_col
		
		# Overwrite existing variable in list
		variable_list[[i]] <- var
	}
	
	marginals <- unique(marginals)
	out <- data.frame()
	for (i in variable_list) {
		temp <- i@output
		missing <- marginals[!marginals %in% colnames(temp)]
		for (j in missing) {
			temp[[j]] <- NA
		}
		out <- rbind(out, temp[c(marginals, new_source_col, new_res_col)])
	}
	
	out <- Ovariable(name, output = out, marginal = c(rep(TRUE, ncol(out) - 2), FALSE, FALSE)) # all marginals except new source and res
	return(out)
}