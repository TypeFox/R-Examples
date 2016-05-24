CombineReadsAcrossRuns <-
function(seqs) {
	## seqs should be a LIST of seq objects with the same format
	## Each element of each seq should be a vector
	## Return: combSeq, a list of the same format as each seqs
	nSeq = length(seqs)
	nVar = length(names(seqs[[1]]))
	combSeq = vector("list", nVar)
	names(combSeq) = names(seqs[[1]])
	for(i in 1:nVar) {
		for(j in 1:nSeq) {
			combSeq[[i]] = c(combSeq[[i]], as.numeric(seqs[[j]][[i]]))
		}
	}
	return(combSeq)
}

