## load library
require("GMD")
require(cluster)

## compute distance using Euclidean metric (default)
data(ruspini)
x <- gdist(ruspini)

## see a dendrogram result by hierarchical clustering
dev.new(width=12, height=6)
plot(hclust(x),
     main="Cluster Dendrogram of Ruspini data",
     xlab="Observations")

## convert to a distance matrix
m <- as.matrix(x)

## convert from a distance matrix
d <- as.dist(m)
stopifnot(d == x)

## Use correlations between variables "as distance"
data(USJudgeRatings)
dd <- gdist(x=USJudgeRatings,method="correlation.of.variables")
dev.new(width=12, height=6)
plot(hclust(dd),
     main="Cluster Dendrogram of USJudgeRatings data",
     xlab="Variables")
