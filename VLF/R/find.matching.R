find.matching <-
function(NucleotideList, AminoAcidList, nuclength, aalength){
	matchAA <- matrix(0, nrow = nrow(AminoAcidList), ncol = aalength+2)
	matchNuc <- matrix(0, nrow = nrow(NucleotideList), ncol = nuclength+2)
	p = 1

	for(i in 1:nrow(NucleotideList)){
		for(n in 1:nrow(AminoAcidList)){
			if(NucleotideList[i,1] == AminoAcidList[n,1]){
				matchAA[p,] <- AminoAcidList[n,]
				matchNuc[p,] <- NucleotideList[i,]
				p = p + 1
			}
		}
	}

	matchAA <- matchAA[1:(p-1),]
	matchNuc <- matchNuc[1:(p-1),]
	
	bothMatch <- list(matchAA, matchNuc)
	
	return(bothMatch)
}
