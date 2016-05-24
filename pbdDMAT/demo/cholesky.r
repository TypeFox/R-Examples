### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

# Initialize process grid
library(pbdDMAT, quiet=T)

if(comm.size() != 2)
  comm.stop("Exactly 2 processors are required for this demo.")

init.grid()

# Setup for the remainder
comm.set.seed(1234, diff=TRUE)
M <- 16
N <- 4
BL <- 2 # blocking --- passing single value BL assumes BLxBL blocking
dA <- ddmatrix("rnorm", M, N, mean=100, sd=10)

A <- as.matrix(dA)

# Cholesky
ch1 <- chol(t(A) %*% A)
ch2 <- as.matrix(chol(t(dA) %*% dA))
comm.print(sum(ch1 - ch2))

# Finish
finalize()
