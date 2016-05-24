### SHELL> mpiexec -np 4 Rscript --vanilla [...].r

# Initial MPI.
library(pbdDEMO, quietly = TRUE)
init.grid()

# Number of rows/columns
n <- 250
p <- 50


# Generate locally only what is needed.
# This will produce a different matrix because of the block cyclic 
# distribution.
comm.set.seed(1234, diff = TRUE)
dx <- ddmatrix.local("rnorm", nrow=n, ncol=p, bldim=4)

print(dx)

finalize()
