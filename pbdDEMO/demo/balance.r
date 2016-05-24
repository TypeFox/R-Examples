### SHELL> mpiexec -np 4 Rscript --vanilla [...].r

### Setup environment.
library(pbdDEMO, quietly = TRUE)
if(comm.size() != 4){
  comm.stop("This example requries 4 processors.")
}
comm.set.seed(1234, diff = TRUE)

### X.gbd can be readed from .csv files distributedly.
N.gbd <- 1 + comm.rank()
X.gbd <- matrix(rnorm(N.gbd * 3), ncol = 3)
comm.cat("X.gbd on rank 2:\n", quiet = TRUE)
comm.print(X.gbd, rank.print = 2, quiet = TRUE)
comm.cat("X.gbd on rank 3:\n", quiet = TRUE)
comm.print(X.gbd, rank.print = 3, quiet = TRUE)

### Run
bal.info <- balance.info(X.gbd)
new.X.gbd <- load.balance(X.gbd, bal.info)
org.X.gbd <- unload.balance(new.X.gbd, bal.info)

comm.cat("\nnew.X.gbd on rank 1:\n", quiet = TRUE)
comm.print(new.X.gbd, rank.print = 1, quiet = TRUE)
comm.cat("\nnew.X.gbd on rank 2:\n", quiet = TRUE)
comm.print(new.X.gbd, rank.print = 2, quiet = TRUE)

comm.cat("\nbal.info on rank 2:\n", quiet = TRUE)
comm.print(bal.info, rank.print = 2, quiet = TRUE)

if(any(org.X.gbd - X.gbd != 0)){
  cat("Unbalance fails in the rank ", comm.rank(), "\n")
}

### Quit.
finalize()
