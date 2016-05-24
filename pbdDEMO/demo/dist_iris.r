# File name: dist_iris.r
# Run: mpiexec -np 4 Rscript dist_iris.r

rm(list = ls())                                       # Clean environment
library(pbdMPI, quietly = TRUE)                       # Load library
library(ape, quietly = TRUE)                          # Load library
if(comm.size() != 4)
  comm.stop("4 processors are required.")

### Load data
X <- as.matrix(iris[, -5])                            # Dimension 150 by 4
X.cid <- iris[, 5]                                    # True id

### Distribute data
jid <- get.jid(nrow(X))
X.gbd <- X[jid,]                                      # GBD row-major format

### Compute distance
X.dist.common <- comm.dist(X.gbd)
X.dist.gbd <- comm.dist(X.gbd, return.type = "gbd")

### Verify
X.dist.gbd <- unlist(allgather(X.dist.gbd[, 3]))
n.diff <- sum(X.dist.common != X.dist.gbd)

### Hierarchical clustering
X.hc <- hclust(X.dist.common, method = "average")
if(comm.rank() == 0){
  X.hc$labels <- as.character(X.cid)
  X.phylo <- as.phylo(X.hc)
  pdf("dist_iris.pdf", width = 10, height = 4)
    plot(X.phylo, tip.color = as.numeric(X.cid) + 1,
         direction = "downwards", no.margin = TRUE,
         label.offset = 0.01, y.lim = c(-0.3, 2.1), cex = 0.5)
  dev.off()

  cat("Difference of common and gbd: ", n.diff, "\n", sep = "")
}

### Finish
finalize()
