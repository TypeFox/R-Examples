slideNucDiag <- 
function(DNAbin, sppVector, width, interval = 1){
	nd <- nucDiag(DNAbin, sppVector)
	if(interval == "codons") interval <- 3
	win <- seq(1, dim(DNAbin)[2] - width, by = interval)
	mat <- matrix(NA, nrow = length(nd), ncol = length(win))
	for(i in 1:length(win)) mat[ ,i] <- sapply(nd, function(x) length(which(x %in% win[i]:(win[i] + width))))
	dimnames(mat)[[1]] <- unique(sppVector)
	mat
}

