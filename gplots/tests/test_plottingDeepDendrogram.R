library(gplots)

dat <- as.matrix(read.csv(file="dat.csv", row.names=1))
dist2 <- function(x) as.dist(1-cor(t(x), method="pearson"))
hclust1 <- function(x) hclust(x, method = "single")

distance <- dist2(dat)
cluster  <- hclust1(distance)

dend <- as.dendrogram(cluster)

## R's default recursion limits will be exceeded when plotting this dendrogram
try( gplots:::plot.dendrogram(dend) )
try( heatmap.2(dat, Rowv=dend) )

## Increase them and try again
options("expressions"=20000)
gplots:::plot.dendrogram(dend)
heatmap.2(dat, Rowv=dend)
