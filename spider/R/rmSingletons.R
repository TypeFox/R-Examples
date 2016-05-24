rmSingletons <- function(sppVector, exclude = TRUE){
	singletons <- names(which(table(sppVector) == 1))
	if(exclude) which(!sppVector %in% singletons) else which(sppVector %in% singletons)
}