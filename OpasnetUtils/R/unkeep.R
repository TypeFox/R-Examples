unkeep <- function (X, cols = NA, sources = FALSE, prevresults = FALSE) 
{
	rescol <- paste(X@name, "Result", sep = "")
	if (sources) 
		cols <- c(cols, colnames(X@output)[grepl("Source$", colnames(X@output))])
	
	for(i in cols) { # Give warning about unkept indices with >1 locations.
		if(
				length(unique(X@output[[i]])) > 1 & 
				X@marginal[match(i, colnames(X@output))]
				) {
			warning(paste(
							"There is >1 unique locations in column", 
							i, ":", 
							paste(as.character(unique(X@output[[i]])), collapse = ", ")
					))
		}
	}
	
	if (prevresults) 
		cols <- c(cols, colnames(X@output)[grepl("Result$", colnames(X@output)) & 
								colnames(X@output) != rescol])
	
	marginals <- colnames(X@output)[X@marginal]
	X@output <- X@output[!colnames(X@output) %in% cols]
	X@marginal <- colnames(X@output) %in% marginals
	
	return(X)
}