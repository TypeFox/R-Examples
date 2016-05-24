### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

# Initialize process grid
library(pbdDMAT, quiet=T)

if(comm.size() != 2)
  comm.stop("Exactly 2 processors are required for this demo.")

init.grid()

# Setup for the remainder
comm.set.seed(25, diff=FALSE)
M <- N <- O <- 16
BL <- 2 # blocking --- passing single value BL assumes BLxBL blocking

A <- matrix(rnorm(M * N, mean = 100, sd = 10), nrow = M, ncol = N)
B <- matrix(rnorm(N * O, mean = 100, sd = 10), nrow = N, ncol = O)

# Distributing matrices
dA <- as.ddmatrix(A, BL)
dB <- as.ddmatrix(B, BL)

# Solve system AX=B and check the serial and parallel results
sol <- solve(A, B)
dsol <- as.matrix(solve(dA, dB))
comm.print(all.equal(sol, dsol))

# Invert matrix A and check the serial and parallel results
inv <- solve(A)
dinv <- as.matrix(solve(dA))
comm.print(all.equal(inv, dinv))

# Finish
finalize()
