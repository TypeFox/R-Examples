mostDE <- function(results, n = 10) {
	if (n < length(results$DEgenes[,1])) {
		topDE = results$DEgenes[1:n,];
	} else {
		sortedGenes = results$allgenes[order(results$allgenes[["p.val.Adjusted"]]), ]
		topDE = sortedGenes[1:n, c(1,3)];
	}
	row.names(topDE) = NULL;
	return(topDE)	
}
