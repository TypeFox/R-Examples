### Setup environment.
suppressMessages(library(pmclust, quietly = TRUE))

### Load data
X <- as.matrix(iris[, -5])

### Standardized
X.std <- scale(X)

### Clustering
suppressMessages(library(pmclust, quietly = TRUE))
comm.set.seed(123, diff = TRUE)

ret.mb1 <- pmclust(X.std, K = 3, method.own.X = "common")
comm.print(ret.mb1)

ret.kms <- pkmeans(X.std, K = 3, method.own.X = "common")
comm.print(ret.kms)

### Finish
finalize()
