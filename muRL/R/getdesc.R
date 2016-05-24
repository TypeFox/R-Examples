getdesc <- function(data){
	
	tmp <- data[(grep("Complete Position Description", data) + 1) : (grep(">Job Listing Duration<", data)-1)]	
	
	tmp <- paste(tmp, collapse = " ")
	return(tmp)
}