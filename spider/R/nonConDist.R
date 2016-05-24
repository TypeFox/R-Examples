nonConDist <-
function(distobj, sppVector = NULL, propZero = FALSE, rmNA = FALSE){
	distobj <- as.matrix(distobj)
	if(length(sppVector) > 0) dimnames(distobj)[[1]] <- sppVector
	nonSpecDists <- list()
	for(i in 1:length(dimnames(distobj)[[1]])){
	  nonSpec <- dimnames(distobj)[[1]] != dimnames(distobj)[[1]][i]
	  nonSpecDists[[i]] <- min(distobj[nonSpec,i] , na.rm = rmNA)
	}	
if(propZero) output <- length(which(unlist(nonSpecDists) == 0))/length(unlist(nonSpecDists)) else output <- unlist(nonSpecDists)

output
}

