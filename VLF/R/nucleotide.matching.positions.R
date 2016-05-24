nucleotide.matching.positions <-
function(matchNuc, nuclength){
	positionNuc <- list()
	p <- 1
	for(i in 1:nrow(matchNuc)){
		for(n in 1:nuclength){
			if(is.na(matchNuc[i,n+2]) == FALSE){
				positionNuc[[p]] <- c(matchNuc[i,1], matchNuc[i, 2], n)
				p = p + 1
			}
		}
	}
	return(positionNuc)
}
