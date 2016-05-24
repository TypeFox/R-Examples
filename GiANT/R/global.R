global.overrepresentation <- function(
		dat,
		geneSet,
		coreSet){
	
	nms <- rownames(dat)

	##########################################
	# 2 by 2 table
	##########################################
	tp <- sum(is.element(geneSet, coreSet))
	# if tp == 0 the p-Value will be 1. No further calculations are necessary
	if(tp == 0){
		return(list(
			pValue = 1,
			coreSet = coreSet,
			intersectGeneSetCoreSet = NULL,
			nAllGenes = length(nms),
			res.all = NULL))
	}
	##########################################
	tn <- sum(!is.element(nms, c(coreSet, geneSet)))
	fp <- length(geneSet) - tp
	fn <- sum(!is.element(coreSet, geneSet))

	table_2by2 <- matrix(c(tp, fp, fn, tn),2,2)

	##########################################
	# R fisher test
	##########################################
	f.res <- fisher.test(table_2by2, alternative = "greater")

	return(list(
		pValue = f.res$p.value,
		coreSet = coreSet,
		intersectGeneSetCoreSet = intersect(geneSet, coreSet),
		nAllGenes = length(nms),
		res.all = f.res))
}

global.ancova <- function(
		dat,
		geneSet,
		labs,
		...){
	requireNamespace("GlobalAncova")

	res <- do.call("GlobalAncova", c(list(xx = dat, test.genes = geneSet, group = labs), ...))

	return(list(
		pValue = res$test.result[2],
		res.all = res))
}

global.test <- function(
		dat,
		geneSet,
		labs,
		...){
	requireNamespace("globaltest")

	res <- do.call("gt", c(list(alternative = t(dat), subsets = geneSet, response = labs), ...))

	return(list(
		pValue = attributes(res)$result[1],
		res.all = res))
}
