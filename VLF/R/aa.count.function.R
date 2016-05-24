aa.count.function <-
function(aminoAcids, seqlength){
	aa <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
	spec.no <- nrow(aminoAcids)
	count <- mat.or.vec(nr = 20, nc = seqlength)
	for(i in 1:seqlength){
		count[1,i] <- length(which(aminoAcids[,i+2] == aa[1]))
		count[2,i] <- length(which(aminoAcids[,i+2] == aa[2]))
		count[3,i] <- length(which(aminoAcids[,i+2] == aa[3]))
		count[4,i] <- length(which(aminoAcids[,i+2] == aa[4]))
		count[5,i] <- length(which(aminoAcids[,i+2] == aa[5]))
		count[6,i] <- length(which(aminoAcids[,i+2] == aa[6]))
		count[7,i] <- length(which(aminoAcids[,i+2] == aa[7]))
		count[8,i] <- length(which(aminoAcids[,i+2] == aa[8]))
		count[9,i] <- length(which(aminoAcids[,i+2] == aa[9]))
		count[10,i] <- length(which(aminoAcids[,i+2] == aa[10]))
		count[11,i] <- length(which(aminoAcids[,i+2] == aa[11]))
		count[12,i] <- length(which(aminoAcids[,i+2] == aa[12]))
		count[13,i] <- length(which(aminoAcids[,i+2] == aa[13]))
		count[14,i] <- length(which(aminoAcids[,i+2] == aa[14]))
		count[15,i] <- length(which(aminoAcids[,i+2] == aa[15]))
		count[16,i] <- length(which(aminoAcids[,i+2] == aa[16]))
		count[17,i] <- length(which(aminoAcids[,i+2] == aa[17]))
		count[18,i] <- length(which(aminoAcids[,i+2] == aa[18]))	
		count[19,i] <- length(which(aminoAcids[,i+2] == aa[19]))
		count[20,i] <- length(which(aminoAcids[,i+2] == aa[20]))
	}	
	rownames(count) <- aa
	return(count)
}
