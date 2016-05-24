SegSeqResProcess <-
function(filename) {
	segseqRaw = read.delim(filename)
	chrs = sort(unique(segseqRaw[,1]))
	segseqRes = vector("list", length(chrs))
	for(j in 1:length(chrs)) {
		segseqRes[[j]] = as.matrix(segseqRaw[segseqRaw[,1]==chrs[j],-1])
	}
	names(segseqRes) = as.character(chrs)
	nChrs = length(segseqRes)
	tauHat = vector("list", nChrs)
	for(j in 1:nChrs) {
		tauHat[[j]] = unique(c(segseqRes[[j]][1,1],segseqRes[[j]][,2]))
	}
	names(tauHat) = names(segseqRes)
	return(tauHat)
}