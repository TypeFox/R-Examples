## Simple dendrogram demo.
##

library(idendr0) # idendro

data(iris)

# perform hierarchical clustering
hc <- hclust(dist(iris[, 1:4]))

# visualize clusters
idendro(hc)
