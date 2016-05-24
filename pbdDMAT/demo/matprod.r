### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

# Initialize process grid
library(pbdDMAT, quiet=T)

if(comm.size() != 2)
  comm.stop("Exactly 2 processors are required for this demo.")

init.grid()

# Generate a random matrix common to all processes and distribute it.
# This approach should only be used while learning the pbdDMAT package.
comm.set.seed(1234, diff=TRUE)
dx <- ddmatrix("rnorm", nrow = 25, ncol = 4)
x <- as.matrix(dx)

# Matrix operations
myprod <- t(dx) %*% dx

# Return the results to a global matrix and compare with R's solution
myprod <- as.matrix(myprod)
rprod <- t(x) %*% x
comm.print(sum(myprod - rprod))

# Finish
finalize()
