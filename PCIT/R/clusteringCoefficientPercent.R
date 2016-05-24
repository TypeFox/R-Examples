clusteringCoefficientPercent <- function(adj) {
	d <- adj[upper.tri(adj)]
	
	d1 <- d[which(d != 0)]
	return(length(d1) / length(d) *100)
}
