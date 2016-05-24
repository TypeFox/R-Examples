### SHELL> mpiexec -np 4 Rscript --vanilla [...].r

### Initial MPI.
library(pbdDEMO, quietly = TRUE)
init()

### Generate balanced fake data.
comm.set.seed(1234, diff = TRUE)
N <- 100                  # Pretend N is large.
p <- 2
### Distributed data.
X.gbd <- matrix(rnorm(N * p), ncol = p)
beta <- 1:p
epsilon <- rnorm(N)
y.gbd <- X.gbd %*% beta + epsilon 

### Run.
ret.gbd <- mpi.ols(y.gbd, X.gbd)

### Output.
comm.print(ret.gbd)
finalize()
