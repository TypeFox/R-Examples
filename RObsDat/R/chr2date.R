# convert a character to date
chr2date <- function(string, tz){
	string <- as.character(string)
	if(all(nchar(string) %in% c(10,19))){
		result <- as.POSIXct(rep(NA, length(string)))
		if(tz!="GMT") warning("Unsupported tz handling other than GMT in chr2date") 
		do.these <- nchar(string) == 10
		result[do.these] <- as.POSIXct(strptime(string[do.these], "%Y-%m-%d", tz="GMT"))
		do.these <- nchar(string) == 19
		result[do.these] <- as.POSIXct(strptime(string[do.these], "%Y-%m-%d %H:%M:%s", tz="GMT"))
		if(any(is.na(result))){
			cat("Unexpected NA in string conversion results\n")
			browser()
		}
	} else {
		cat("Unimplemented date format in chr2date")
		browser()
	}
	if(any(is.null(result))){
		cat("Invalid Date conversion in chr2date. Entering browser mode\n")
		browser()
	}
	return(result)
}
