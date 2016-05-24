# File name: iris_gbd.r
# Run: mpiexec -np 4 Rscript iris_gbd.r

rm(list = ls())                                       # Clean environment
library(pbdMPI, quietly = TRUE)                         # Load library
if(comm.size() != 4)
  comm.stop("4 processors are required.")

### Load data
X <- as.matrix(iris[, -5])                            # Dimension 150 by 4
X.cid <- as.numeric(iris[, 5])                        # True id

### Distribute data
jid <- get.jid(nrow(X))
X.gbd <- X[jid,]                                      # GBD row-major format

### Standardized
N <- allreduce(nrow(X.gbd))                           # 150
p <- ncol(X.gbd)                                      # 4
mu <- allreduce(colSums(X.gbd / N))
X.std <- sweep(X.gbd, 2, mu, FUN = "-")               # Substract mean
std <- sqrt(allreduce(colSums(X.std^2 / (N - 1))))
X.std <- sweep(X.std, 2, std, FUN = "/")              # Divide standard error

### SVD manually in serial
X.tmp <- crossprod(X.std)                             # X'X (local)
X.tmp <- allreduce(X.tmp)
dim(X.tmp) <- c(p, p)
ret <- eigen(X.tmp)                                   # X'X = V D^2 V'
d <- sqrt(ret$values)
v <- ret$vectors
u <- X.std %*% v %*% diag(1/d)                        # Why X V D^(-1)) = U?

### Validation 
tmp <- svd(scale(X))
comm.print(c(sum(abs(tmp$u[jid,] - u)),
             sum(abs(tmp$v - v)),
             sum(abs(tmp$d - d))))

### Project on column space of singular vectors
A <- u %*% diag(d)
B <- X.std %*% v                                      # A ~ B 
X.prj <- A[, 1:2]                                     # Only useful for plot
X.prj <- do.call("rbind", allgather(X.prj))

### Clustering
library(pmclust, quietly = TRUE)
comm.set.seed(123, diff = TRUE)

X.spmd <- X.std
PARAM.org <- set.global(K = 3)                        # Preset storage
.pmclustEnv$CONTROL$debug <- 0                        # Disable debug messages
PARAM.org <- initial.center(PARAM.org)                # Initial parameters
PARAM.kms <- kmeans.step(PARAM.org)                   # K-means
X.kms.cid <- allgather(.pmclustEnv$CLASS.spmd,
                       unlist = TRUE)

PARAM.org <- set.global(K = 3)                        # Preset storage
.pmclustEnv$CONTROL$debug <- 0                        # Disable debug messages
PARAM.org <- initial.em(PARAM.org,
                        MU = PARAM.kms$MU)            # Initial by K-means
PARAM.mbc1 <- em.step(PARAM.org)                      # Model-based clustering
X.mbc1.cid <- allgather(.pmclustEnv$CLASS.spmd,
                        unlist = TRUE)

PARAM.org <- set.global(K = 3, RndEM.iter = 1000)     # Preset storage
.pmclustEnv$CONTROL$debug <- 0                        # Disable debug messages
PARAM.org <- initial.RndEM(PARAM.org)                 # Initial by Rand-EM
PARAM.mbc2 <- em.step(PARAM.org)                      # Model-based clustering
X.mbc2.cid <- allgather(.pmclustEnv$CLASS.spmd,
                        unlist = TRUE)

### Validation
X.kms.adjR <- EMCluster::RRand(X.cid, X.kms.cid)$adjRand
X.mbc1.adjR <- EMCluster::RRand(X.cid, X.mbc1.cid)$adjRand
X.mbc2.adjR <- EMCluster::RRand(X.cid, X.mbc2.cid)$adjRand
comm.print(c(X.kms.adjR, X.mbc1.adjR, X.mbc2.adjR))

### Swap classification id
tmp <- X.kms.cid
X.kms.cid[tmp == 1] <- 2
X.kms.cid[tmp == 2] <- 1
tmp <- X.mbc1.cid
X.mbc1.cid[tmp == 1] <- 2
X.mbc1.cid[tmp == 2] <- 1

### Display on first 2 components
if(comm.rank() == 0){
  pdf("gbd_plot.pdf")
  
  par(mfrow = c(2, 2))
  plot(X.prj, col = X.cid + 1, pch = X.cid,
       main = "iris (true)", xlab = "PC1", ylab = "PC2")
  plot(X.prj, col = X.kms.cid + 1, pch = X.kms.cid,
       main = paste("iris (kmeans)", sprintf("%.4f", X.kms.adjR)),
       xlab = "PC1", ylab = "PC2")
  plot(X.prj, col = X.mbc1.cid + 1, pch = X.mbc1.cid,
       main = paste("iris (model-based 1)", sprintf("%.4f", X.mbc1.adjR)),
       xlab = "PC1", ylab = "PC2")
  plot(X.prj, col = X.mbc2.cid + 1, pch = X.mbc2.cid,
       main = paste("iris (model-based 2)", sprintf("%.4f", X.mbc2.adjR)),
       xlab = "PC1", ylab = "PC2")
  
  dev.off()
}

### Finish
finalize()
