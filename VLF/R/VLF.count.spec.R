VLF.count.spec <-
function (freq, p,seqlength){
	count <- mat.or.vec(nr = nrow(freq), nc = 1)
	for (n in 1:nrow(freq)){
		count[n] <- length(which(freq[n,] < p))
	}
	return(count)
}
