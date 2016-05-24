### SHELL> mpiexec -np 4 Rscript --vanilla [...].r

### Initial MPI.
library(pbdDEMO, quietly = TRUE)
init()

### Generate balanced fake data.
comm.set.seed(1234)
N <- 100 * comm.size()    # Pretend N is large.
x <- rnorm(N)
id.get <- get.jid(N)
x.gbd <- x[id.get]        # Distributed data.

### Run.
ret.gbd <- mpi.stat(x.gbd)

### Output.
comm.print(ret.gbd)
finalize()
