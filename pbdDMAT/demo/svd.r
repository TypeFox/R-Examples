### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

# Initialize process grid
library(pbdDMAT, quiet=T)

if(comm.size() != 2)
  comm.stop("Exactly 2 processors are required for this demo.")

init.grid()

# Setup for the remainder
comm.set.seed(1234, diff=TRUE)
M <- N <- 16
BL <- 2 # blocking --- passing single value BL assumes BLxBL blocking

dA <- ddmatrix("rnorm", nrow=M, ncol=N, mean=100, sd=10)
A <- as.matrix(dA)

# SVD
svd1 <- La.svd(A)
svd2 <- La.svd(dA)
svd2 <- lapply(svd2, as.matrix)

# Test equality of serial and parallel versions
comm.print(all.equal(svd1$d, as.vector(svd2$d)))
comm.print(all.equal(svd1$u, svd2$u))
comm.print(all.equal(svd1$vt, svd2$vt))

# Finish
finalize()
