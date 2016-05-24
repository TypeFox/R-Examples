search.BOLD <- 
function(taxon, exhaustive = FALSE){
	space <- grep(" ", taxon)
	if(length(space) > 0) taxon[space] <- sapply(space, function(x) paste("\"", taxon[x], "\"", sep=""))
	taxon <- paste(taxon, collapse="%20")
	taxon <- gsub(" ", "%20", taxon)
	if(exhaustive){
		num <- stats.BOLD(taxon)
		if(num > 500){
			offsetVals <- 500 * 0:(as.integer(num/500))
			IDs <- list()
			for(i in 1:length(offsetVals)){
				URL <- paste("http://v3.boldsystems.org/index.php/Public_Ajax_RecordList?query=", taxon, "&limit=500&offset=", offsetVals[i], sep="")
				res <- scan(file = URL, what = "", sep = "\n", quiet = TRUE)
				res <- res[grep("selected_processids", res)]
				IDs[[i]] <- sapply(strsplit(res, "\""), function(x) x[6])
			}
			IDs <- unlist(IDs)
			} else exhaustive <- FALSE
		}
	if(!exhaustive){
		URL <- paste("http://v3.boldsystems.org/index.php/Public_Ajax_RecordList?query=", taxon, "&limit=500&offset=0", sep="")
		res <- scan(file = URL, what = "", sep = "\n", quiet = TRUE)
		res <- res[grep("selected_processids", res)]
		IDs <- sapply(strsplit(res, "\""), function(x) x[6])
		}
	if(length(IDs) == 0) stop("No taxa with that name found on BOLD") 
		else if(length(IDs) == 500){
			warning("Maximum number of records reached. Try using 'exhaustive = TRUE'")
			IDs
		} else IDs
}