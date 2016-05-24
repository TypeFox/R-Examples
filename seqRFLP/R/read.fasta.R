read.fasta <-
function(file = NULL){
    if(is.null(file))
	{stop("Please specify the input fasta file.")}
	result <- readLines(file)
	nameline <- result[grepl("[>]", result)]
	test <- regexpr(">", nameline) > 1
	if(any(test)){
	warning(paste("\">\" in line(s)",which(test),"\n appeared not at the beginning. \n Please remove any character(s) before \">\"."))
	}
	result = result[grepl("[A-Za-z0-9]", result)]
	result <- ConvFas(result, "fas")
	if(any(regexpr(">", result[seq(1, length(result), by = 2)]) < 0)){
	    xx <-  2 * which(regexpr(">", result[seq(1, length(result), by = 2)]) < 0)
		if( length(xx) > 10){ 
		   xx <- xx[1:10]
		   stop(paste("readfasta could not find \">\" in row: \n",
			    paste(xx,  collapse = ", "),"... \n",
	            "Make sure the file ", file, " is in fasta format.\n"))
		   }
	   }
	class(result) <- "fasta"
	return(result)
}

