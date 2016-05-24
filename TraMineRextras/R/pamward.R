## runs a pam from the solution in k groups of a hierarchical clustering
## author: Gilbert Ritschard

pamward <- function(dist, k=3, method="ward") {
	clustw <- agnes(dist, diss=T, method=method)
	clw <- cutree(clustw, k=k)
	centers <- disscenter(dist, group=clw, medoids.index="first")
	return <- pam(dist, diss=T, k=k, medoids=centers)
}

