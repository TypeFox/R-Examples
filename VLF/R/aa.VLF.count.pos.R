aa.VLF.count.pos <-
function(freq, p, seqlength){
	count <- mat.or.vec(nr = 1, nc = seqlength)
	for(i in 1:seqlength){
		count[i] <- length(which(freq[,i] < p))
	}
	return(count)
}
