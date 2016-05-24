library(rbenchmark)
library(fastcluster)
library(Rclusterpp)

ROWS    <- c(50, 100, 500, 1000, 2000)
COLUMNS <- 10;

results <- c()

for (r in ROWS) {
	data <- matrix(rnorm(r * COLUMNS), nrow=r)
	result <- benchmark( 
		Rclusterpp = Rclusterpp.hclust(data, method="ward"),
		hclust = stats::hclust((dist(data, method="euclidean")^2)/2.0, method="ward"),
		fastcluster = fastcluster::hclust((dist(data, method="euclidean")^2)/2.0, method="ward"),
		replications = 5, 
		columns=c("test", "elapsed", "user.self", "sys.self"),
		order="elapsed"
	)
	results <- rbind(results, cbind(result, obs = rep(r, nrow(result)), method = rep("ward", nrow(result))))
}
print(results)

results <- c()
for (r in ROWS) {
	data <- matrix(rnorm(r * COLUMNS), nrow=r)
	result <- benchmark( 
		Rclusterpp = Rclusterpp.hclust(data, method="average", distance="euclidean"),
		RclusterppDistance = Rclusterpp.hclust(dist(data, method="euclidean"), method="average"), 
		hclust = stats::hclust(dist(data, method="euclidean"), method="average"),
		fastcluster = fastcluster::hclust(dist(data, method="euclidean"), method="average"),
		replications = 5, 
		columns=c("test", "elapsed", "user.self", "sys.self"),
		order="elapsed"
	)
	results <- rbind(results, cbind(result, obs = rep(r, nrow(result)), method = rep("average:euclidean", nrow(result))))
}
print(results)
