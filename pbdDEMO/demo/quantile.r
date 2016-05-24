### SHELL> mpiexec -np 4 Rscript --vanilla [...].r

### Initial MPI.
library(pbdDEMO, quietly = TRUE)
init()

### Generate balanced fake data.
comm.set.seed(1234, diff = TRUE)
N <- 100                  # Pretend N is large.
x.gbd <- rnorm(N)         # Distributed data.

### Run.
ret.gbd <- mpi.quantile(x.gbd)

### Output.
comm.print(ret.gbd)
finalize()
