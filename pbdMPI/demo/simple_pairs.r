### SHELL> mpiexec -np 4 Rscript --vanilla [...].r

### Initial.
suppressMessages(library(pbdMPI, quietly = TRUE))

### Original distance.
X <- matrix(unlist(iris[, -5]), ncol = 4)
X.dist.org <- dist(X)

### Get a gbd row-block matrix.
# if(comm.rank() != 0){
#   X <- matrix(0, nrow = 0, ncol = 4)
# }
X.gbd <- comm.as.gbd(X)

### Get a common distance matrix.
X.dist.common <- comm.dist(X.gbd)

### Get a distributed distance object.
X <- matrix(unlist(iris[, -5]), ncol = 4)
pairsid.gbd <- comm.allpairs(nrow(X))
X.dist.gbd.1 <- comm.pairwise(X, pairsid.gbd)

### The other way.
X.dist.gbd.2 <- comm.pairwise(X.gbd)

### Verify.
d.org <- as.vector(X.dist.org)
d.1 <- as.vector(X.dist.common)
d.2 <- do.call("c", allgather(X.dist.gbd.1[, 3]))
d.3 <- do.call("c", allgather(X.dist.gbd.2[, 3]))
comm.print(all(d.org == d.1))
comm.print(all(d.org == d.2))
comm.print(all(d.org == d.3))

### Finish.
finalize()
