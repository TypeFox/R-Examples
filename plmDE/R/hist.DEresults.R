hist.DEresults <- function(x, ...) {
	pValDistribution = x$allgenes[["p.val.Adjusted"]]
	hist(pValDistribution, main = NULL, xlab = "Adjusted p-values", freq = FALSE, col = colors()[311], ...); box()
}
