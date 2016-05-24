## Simple demo showing how to make `idendro' visualize long dimnames
## and case names in the heat map.
##

library(idendr0) # idendro

# generate data in feature space
n <- 10
x <- data.frame(very.long.dim.name1 = c(rnorm(n, -1), rnorm(n, 1)),
    long.dim.name2 = c(rnorm(n, -1), rnorm(n, 1)))
rownames(x) <- paste('my brand new observation', 1:(2*n))

# compute pairwise distances
dx <- dist(x)

# perform hierarchical clustering
hx <- hclust(dx)

# visualize clusters
# setting margins large enough to hold the long dimnames
idendro(hx, x, mai = par('mai') * c(1, 1, 1.7, 4))
