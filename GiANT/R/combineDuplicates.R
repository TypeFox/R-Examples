##########################################################################
#mergeProbesForGenes
#
# dat:		Matrix of data with one row per gene.
# method: 	The combination method for the rows related to the same gene.
#
# Computes a combination of duplicated entries in an expression matrix.
# A duplication free expression matrix is returned.
##########################################################################
mergeProbesForGenes <- function(dat, method = c("mean", "max", "min", "median")){
	
	method <- match.arg(method)

	method <- getFunction(method)
	nms <- unique(rownames(dat))

	tmp <- t(sapply(nms, function(x){apply(dat[rownames(dat) == x,,drop=FALSE], 2, method)}))
	rownames(tmp) <- nms
	colnames(tmp) <- colnames(dat)

	return(tmp)
}