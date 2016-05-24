seqStat <- function(DNAbin, thresh = 500){
	rr <- sapply(DNAbin, FUN = function(x) length(x) - checkDNA(x, gapsAsMissing=TRUE))
	tab <- table(NULL)
	tab[1:5] <- c(min(rr), max(rr), mean(rr), median(rr), length(which(rr < thresh)))
	names(tab) <- c("Min", "Max", "Mean", "Median", "Thresh")
	round(tab, digits=0)
}
