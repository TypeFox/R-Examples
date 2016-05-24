library(QZ, quiet = TRUE)

### Generate Data
set.seed(123)
X <- matrix(rnorm(500), nrow = 25)
X <- t(X) %*% X
A <- X[1:8, 9:20]
B <- X[1:8, 1:8]
C <- X[9:20, 9:20]

### Perform generalized eigenanalysis 
ret.qz <- fda.geigen(A, B, C)
ret.fda <- fda::geigen(A, B, C)

### Verify
round(abs(ret.qz$values - ret.fda$values))
round(abs(ret.qz$Lmat - ret.fda$Lmat))
round(abs(ret.qz$Mmat - ret.fda$Mmat))

