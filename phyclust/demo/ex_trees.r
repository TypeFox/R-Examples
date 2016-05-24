library(phyclust, quiet = TRUE)

### Examples to draw Neighbor-Joining trees.
set.seed(1234)
K <- 3                      # Number of clusters
N <- 50                     # Number of sequences

N.K <- c(15, 20, 15)        # Number of sequences in clusters.
X.class <- rep(1:K, N.K)    # Cluster id for each sequence.

### Generate trees in different rates.
tree.K.1 <- gen.unit.K(K, N.K, rate.anc = 1, rate.dec = 10)
tree.K.2 <- gen.unit.K(K, N.K, rate.anc = 10, rate.dec = 1)

### Plot the trees by plot.nj() with type "p".
par(mfrow = c(2, 3))
plotnj(tree.K.1$max, X.class = X.class, type = "p", main = "max")
axis(1)
plotnj(tree.K.1$equal, X.class = X.class, type = "p", main = "equal")
axis(1)
plotnj(tree.K.1$star, X.class = X.class, type = "p", main = "star")
axis(1)
plotnj(tree.K.2$max, X.class = X.class, type = "p", main = "max")
axis(1)
plotnj(tree.K.2$equal, X.class = X.class, type = "p", main = "equal")
axis(1)
plotnj(tree.K.2$star, X.class = X.class, type = "p", main = "star")
axis(1)

