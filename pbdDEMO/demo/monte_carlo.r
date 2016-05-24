### SHELL> mpiexec -np 4 Rscript --vanilla [...].r

### Setup environment.
library(pbdDEMO, quietly = TRUE)
comm.set.seed(1234, diff = TRUE)

### Run
N.gbd <- 1000
X.gbd <- matrix(runif(N.gbd * 2), ncol = 2)
r.gbd <- sum(sqrt(rowSums(X.gbd^2)) <= 1)
ret <- allreduce(c(N.gbd, r.gbd), op = "sum")
PI <- 4 * ret[2] / ret[1]
comm.print(PI)

### Quit.
finalize()
