# File name: iris_dmat.r
# Run: mpiexec -np 4 Rscript iris_dmat.r

rm(list = ls())                                       # Clean environment
library(EMCluster)

### Load data
X <- as.matrix(iris[, -5])                            # Dimension 150 by 4
X.cid <- as.numeric(iris[, 5])                        # True id

### Convert to matrix
X.spmd <- as.matrix(X)

### Standardized
X.std <- scale(X.spmd)

init.EM(X.std, 3)
