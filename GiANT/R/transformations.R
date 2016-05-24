##########################################################################
#transformations
#
# x:			One value per gene gained from gene level statistic.
#
# A collection of methods which can be used to transform the gene level
# statistic.
##########################################################################
transformation.abs <- function(x){
	tmp <- abs(x)
	names(tmp) <- names(x)
	return(tmp)
}

transformation.square <- function(x){
	tmp <- x^2
	names(tmp) <- names(x)
	return(tmp)
}

transformation.localFdr <- function(x,
		statistic="pvalue",
		cutoff.method="fndr",
		pct0=0.75){

	requireNamespace("fdrtool")
	tmp <- fdrtool::fdrtool(x,
		statistic = statistic,
		cutoff.method = cutoff.method,
		pct0 = pct0,
		plot=FALSE,
		verbose = FALSE)$lfdr
	names(tmp) <- names(x)
	return(tmp)
}

transformation.binarize <- function(x, quant = 0.95){
	tmp <- rep(0, length(x))
	tmp[x > quantile(x, prob = quant)] <- 1
	names(tmp) <- names(x)
	return(tmp)	
}

transformation.rank <- function(x){
	tmp <- rank(x)
	names(tmp) <- names(x)
	return(tmp)
}

transformation.adjust <- function(x, adjMethod = "fdr"){
	tmp <- p.adjust(x, method = adjMethod)
	names(tmp) <- names(x)
	return(tmp)
}

transformation.adjustAndBinarize <- function(x, adjMethod = "fdr", threshold = 0.05){
	tmp <- as.numeric(p.adjust(x, method = adjMethod) < threshold)
	names(tmp) <- names(x)
	return(tmp)
}