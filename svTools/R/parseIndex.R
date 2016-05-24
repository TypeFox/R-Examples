parseIndex <- function (file)
{  
	empty <- data.frame(index = character(), description = character(),
		stringsAsFactors = FALSE)
	if (!file.exists(file)) return(empty)
	rl <- readLines(file) 
	if (!length(rl)) return(empty)
	lines <- (regexpr("^", rl) > 0)
	index <- gsub(" +.*$", "", rl[!lines])
	description <- gsub("^.*? +", "", rl[!lines])
	return(data.frame(index = index, description = description,
		stringsAsFactors = FALSE))
}
