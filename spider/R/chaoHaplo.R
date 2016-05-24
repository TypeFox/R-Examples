chaoHaplo <- function(DNAbin){
	haplo <- haplotype(DNAbin)
	i <- if(length(grep("[-|?|r|y|m|k|w|s|b|d|h|v|n]", DNAbin)) > 0) message("There are missing or ambiguous data, which may cause an overestimation of the number of haplotypes")
	nums <- sapply(attr(haplo, "index"), length)
	n <- dim(DNAbin)[1]
	h <- length(nums)
	s <- length(which(nums == 1))
	d <- length(which(nums == 2))
	#Estimated number of haplotypes (From Vink et al 2011)
	if(d > 0) est <- h + ((s^2)/(2 * d)) else est <- h + ((s * (s - 1))/2)
	
	#Confidence intervals (Modified from Chao 1989)
	varest <- est/((h / (est - h)) - n/est)
	C <- exp(1.96 * sqrt(log( 1 + (varest / ((est - h)^2)))))
	
	low <- h + (est - h)/C
	high <- h + (est - h) * C
    
	c(est, low, high)
}
	

