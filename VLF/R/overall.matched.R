overall.matched <-
function(positionNuc, positionAA, nuclength, aalength){
	codonPos <- rep(1:3, nuclength)
	codon <- rep(1:aalength, each = 3)
	
	matched = list()
	p <- 1
	for(i in 1:length(positionAA)){
		for(n in 1:length(positionNuc)){
			if(positionAA[[i]][1] == positionNuc[[n]][1]){
				if(as.numeric(positionAA[[i]][3]) == codon[as.numeric(positionNuc[[n]][3])]){
					matched[[p]] <- c(positionAA[[i]][1], positionAA[[i]][2], positionAA[[i]][4], codon[as.numeric(positionNuc[[n]][3])], codonPos[as.numeric(positionNuc[[n]][3])])
					p = p + 1
				}
			}
		}
	}
	return(matched)
}
