### SHELL> mpiexec -np 4 Rscript --vanilla [...].r

### Initial MPI.
library(pbdDEMO, quietly = TRUE)
init.grid()
if(comm.size() != 4){
  comm.stop("This example requries 4 processors.")
}
comm.set.seed(1234, diff = TRUE)

### X.gbd can be readed from .csv files distributedly.
N.gbd <- 1 + comm.rank()
X.gbd <- matrix(rnorm(N.gbd * 3), ncol = 3)

### Run.
X.dmat <- gbd2dmat(X.gbd)
X <- as.matrix(X.dmat)
new.X.gbd <- dmat2gbd(X.dmat)

### Output.
if(comm.rank() == 1){
  cat("(local,part) new.X.gbd on rank = 1:\n")
  print(new.X.gbd)
}
if(comm.rank() == 2){
  cat("\n(global,all) X[4:6,] on all processors:\n")
  print(X[4:6,])
}
finalize()
