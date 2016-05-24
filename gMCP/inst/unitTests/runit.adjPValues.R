checkAdjPValues <- function(graph, pvalues) {	
	adjP <- gMCP:::adjPValues(graph, pvalues)@adjPValues
	for (p in adjP) {
		result <- gMCP(graph, pvalues, alpha=p)
		checkEquals(adjP<=p,getRejected(result)) 
	}	
}

test.adjPValues <- function() {
	# Bonferroni-Holm
	graph <- BonferroniHolm(3)
	pvalues <- c(0.01, 0.02, 0.05)
	checkAdjPValues(graph, pvalues)
	# Parallel Gatekeeping
	graph <- parallelGatekeeping()
	pvalues <- c(0.01, 0.01, 0.03, 0.04)
	checkAdjPValues(graph, pvalues)
}