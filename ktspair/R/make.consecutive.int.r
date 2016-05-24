make.consecutive.int <- function(y) {
	intialwarn = getOption("warn")
	options(warn = -1)
	if(is.null(y)) {return(NULL)}  
	if(!is.vector(y)){
		y = as.vector(as.character(y))
	}
	out <- as.integer(as.factor(as.character(y)))-1  
	options(warn = intialwarn)
	return(out)
}

