rm(list = ls())                                       # Clean environment


### Load data
X <- as.matrix(iris[, -5])                            # Dimension 150 by 4
X.cid <- as.numeric(iris[, 5])                        # True id


### Transformation and check
X.std <- scale(X)                                     # Standardize
mu <- colMeans(X.std)                                 # Columns means are near 0
cov <- cov(X.std)                                     # Diagonals are near 1
print(mu)
print(cov)


### SVD
X.svd <- svd(X.std)


### Project on column space of singular vectors
A <- X.std %*% diag(X.svd$d)
B <- X.std %*% X.svd$v
C <- prcomp(X.std)$x                                # A = B = C

X.prj <- C[, 1:2]                                   # project onto first 2 PC's


### Clustering
set.seed(1234)                                        # Set overall seed
X.kms <- kmeans(X.std, 3)                             # K-means
X.kms
X.kms.cid <- X.kms$cluster                            # Classification

library(EMCluster)                                    # Model-based clustering
X.mbc <- init.EM(X.std, 3)                            # Initial by em-EM
X.mbc
X.mbc.cid <- X.mbc$class                              # Classification


### Validation
X.kms.adjR <- RRand(X.cid, X.kms.cid)$adjRand       # Adjusted Rand index
X.mbc.adjR <- RRand(X.cid, X.mbc.cid)$adjRand


### Swap classification id
X.kms.cid[X.kms.cid == 2] <- 4
X.kms.cid[X.kms.cid == 3] <- 2
X.kms.cid[X.kms.cid == 4] <- 3


### Display on first 2 components
pdf("serial_plot.pdf")

par(mfrow = c(2, 2))
plot(X.prj, col = X.cid + 1, pch = X.cid,
     main = "iris (true)", xlab = "PC1", ylab = "PC2")
plot(X.prj, col = X.kms.cid + 1, pch = X.kms.cid,
     main = paste("iris (k-Means)", sprintf("%.4f", X.kms.adjR)),
     xlab = "PC1", ylab = "PC2")
plot(X.prj, col = X.mbc.cid + 1, pch = X.mbc.cid,
     main = paste("iris (Model-based)", sprintf("%.4f", X.mbc.adjR)),
     xlab = "PC1", ylab = "PC2")
accuracy <- c(X.kms.adjR, X.mbc.adjR)
names(accuracy) <- c("k-Means", "Model-based")
barplot(accuracy, main = "Clustering Accuracy")

dev.off()
