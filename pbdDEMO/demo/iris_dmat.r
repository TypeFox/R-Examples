# File name: iris_dmat.r
# Run: mpiexec -np 4 Rscript iris_dmat.r

rm(list = ls())                                       # Clean environment
library(pbdDMAT, quietly = TRUE)                        # Load library
init.grid()
#if(comm.size() != 4)
#  comm.stop("4 processors are required.")



dmat_opts$BLDIM <- 5



### Load data
X <- as.matrix(iris[, -5])                            # Dimension 150 by 4
X.cid <- as.numeric(iris[, 5])                        # True id

### Convert to ddmatrix
X.dmat <- as.ddmatrix(X)

### Standardized
X.std <- scale(X.dmat)
#mu <- as.matrix(colMeans(X.std))
#cov <- as.matrix(cov(X.std))
#comm.print(mu)
#comm.print(cov)

### SVD
#X.svd <- svd(X.std)

### Project on column space of singular vectors
#A <- X.svd$u %*% diag(X.svd$d, type="ddmatrix")
#B <- X.std %*% X.svd$v                                # A ~ B
#X.prj <- as.matrix(A[, 1:2])                          # Only useful for plot

### Clustering
library(pmclust, quiet = TRUE)
comm.set.seed(123, diff = TRUE)

X.dmat <- X.std
PARAM.org <- set.global.dmat(K = 3)                   # Preset storage
.pmclustEnv$CONTROL$debug <- 0                        # Disable debug messages
PARAM.org <- initial.center.dmat(PARAM.org)
PARAM.kms <- kmeans.step.dmat(PARAM.org)              # K-means
X.kms.cid <- as.vector(.pmclustEnv$CLASS)

### Validation
X.kms.adjR <- EMCluster::RRand(X.cid, X.kms.cid)$adjRand
comm.print(X.kms.adjR)

### Swap classification id
tmp <- X.kms.cid
X.kms.cid[tmp == 1] <- 3
X.kms.cid[tmp == 2] <- 1
X.kms.cid[tmp == 3] <- 2

### Display on first 2 components
#if(comm.rank() == 0){
#  pdf("dmat_plot.pdf")
#  
#  par(mfrow = c(2, 2))
#  plot(X.prj, col = X.cid + 1, pch = X.cid,
#       main = "iris (true)", xlab = "PC1", ylab = "PC2")
#  plot(X.prj, col = X.kms.cid + 1, pch = X.kms.cid,
#       main = paste("iris (kmeans)", sprintf("%.4f", X.kms.adjR)),
#       xlab = "PC1", ylab = "PC2")
#  dev.off()
#}

### Finish
finalize()

