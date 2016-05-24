library(MixSim, quietly=TRUE)

### iris data
X <- as.matrix(iris[, -5])                                 # Iris matrix
X.cid <- as.numeric(iris[, 5])                             # True ids
N <- nrow(X)                                               # 150 flowers
p <- ncol(X)                                               # 4 dimensions
K <- length(unique(X.cid))                                 # 3 classes

### Empirical parameters based on true ids.
ETA <- as.vector(table(X.cid)) / N                         # Mixing proportion
MU <- rbind(colMeans(X[X.cid == 1,]),                      # Centers
            colMeans(X[X.cid == 2,]),
            colMeans(X[X.cid == 3,]))
S <- c(cov(X[X.cid == 1,]),
       cov(X[X.cid == 2,]),
       cov(X[X.cid == 3,]))
dim(S) <- c(p, p, K)

### Compute overlaps.
(ret <- overlap(ETA, MU, S))
(levels(iris[, 5]))
