# File name: dist_pony.r
# Run: mpiexec -np 4 Rscript dist_pony.r

rm(list = ls())                                       # Clean environment
library(pbdMPI, quietly = TRUE)                       # Load library
library(phyclust, quietly = TRUE)                     # Load library
if(comm.size() != 4)
  comm.stop("4 processors are required.")

### Load data
data.path <- paste(.libPaths()[1], "/phyclust/data/pony524.phy", sep = "")
pony.524 <- read.phylip(data.path)
X <- pony.524$org

### Distribute data
jid <- get.jid(nrow(X))
X.gbd <- X[jid,]                                      # GBD row-major format

### Compute distance
X.dist.common <- phyclust.edist(X)

### User-defined distance
my.dist <- function(x, y){
  as.vector(phyclust.edist(rbind(x, y)))
}
X.gbd.dist <- comm.pairwise(X.gbd, FUN = my.dist)

### Verify
X.dist <- unlist(allgather(X.gbd.dist[, 3]))
n.diff <- sum(X.dist.common != X.dist)

### Neighbor Joining
X.nj <- nj(X.dist.common)
if(comm.rank() == 0){
  set.seed(1234)
  pony.phyclust <- phyclust(X, 3)

  pdf("dist_pony.pdf", width = 6, height = 5)
    par(mar = c(0, 0, 2, 0))
    plotnj(X.nj, X.class = pony.phyclust$class.id)
    title("Pony 524 & K = 3")
  dev.off()

  cat("Difference of common and gbd: ", n.diff, "\n", sep = "")
}

### Finish
finalize()
