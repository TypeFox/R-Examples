clusteringCoefficient <- function(adj, FUN='localClusteringCoefficient', ...) {
	FUN <- match.fun(FUN)
	
	coef <- FUN(adj, ...)
	
	return(coef)
}
