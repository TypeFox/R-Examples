c.skatCohort <- function(...){
  ev <- environment()
  cl <- match.call()
	ncohort <- length(cl)-1
	re <- NULL
	for(i in 1:ncohort){
		re <- c(re,as.list(eval(cl[[i+1]],ev) ))
	}
	class(re) <- "skatCohort"
	re
	}
	
c.seqMeta <- function(...){
	ev <- environment()
  cl <- match.call()
	ncohort <- length(cl)-1
	re <- NULL
	for(i in 1:ncohort){
		re <- c(re,as.list(eval(cl[[i+1]],ev) ))
	}
	class(re) <- "seqMeta"
	re
	}