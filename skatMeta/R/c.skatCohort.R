c.skatCohort <- function(...){
	cl <- match.call()
	ncohort <- length(cl)-1
	re <- NULL
	for(i in 1:ncohort){
		re <- c(re,as.list(eval(cl[[i+1]]) ))
	}
	class(re) <- "skatCohort"
	re
	}