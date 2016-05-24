geturl <- function(callmatch){
	if(!"printurl" %in% names(callmatch) || !isTRUE(callmatch$printurl)){
		return();
	}
	callmatch <- as.list(callmatch);
	FUN <- callmatch[[1]];
	args <- callmatch[-1];
	
	url = "http://rdev.mobilizingcs.org/R/png/Mobilize/";
	url <- paste(url,FUN,"?",sep="");
	
	skipargs <- c("printurl", "verbose")
	
	for(i in 1:length(args)){
		if(names(args[i]) %in% skipargs){
			next;
		}
		url = paste(url, "&", names(args[i]), "='", eval(args[[i]]), "'", sep="")
	}
	print(url);
	return(url);
}