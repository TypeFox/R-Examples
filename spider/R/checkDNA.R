checkDNA <-
function(DNAbin, gapsAsMissing = TRUE){
	if(gapsAsMissing) bases <- c(2, 240, 4) else bases <- c(2, 240)
	if(is.list(DNAbin)) output <- sapply(DNAbin, function(x) length(which(as.numeric(x) %in% bases)))
	if(is.matrix(DNAbin)) output <- apply(DNAbin, MARGIN = 1, FUN = function(x) length(which(as.numeric(x) %in% bases)))
output
}

