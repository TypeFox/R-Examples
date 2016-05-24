aminoAcid.matching.positions <-
function(matchAA, aalength){
	positionAA <- list()
	p <- 1
	for(i in 1:nrow(matchAA)){
		for(n in 1:aalength){
			if(is.na(matchAA[i, n+2]) == FALSE){
				positionAA[[p]] <- c(matchAA[i,1], matchAA[i, 2], n, matchAA[i,n+2])
				p = p + 1
			}
		}
	}
	return(positionAA)
}
