### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

### Initial.
suppressMessages(library(pbdMPI, quietly = TRUE))
init()
.comm.size <- comm.size()
.comm.rank <- comm.rank()

if(.comm.size != 2){
  comm.stop("2 processors are requried.\n", quiet = TRUE)
  finalize()
}

### Examples.
if(.comm.rank == 0){
  a <- c(T, F, NA)
} else{
  a <- T
}

tmp <- comm.any(a)
comm.print(tmp)
tmp <- comm.all(a)
comm.print(tmp)
tmp <- comm.any(a, na.rm = TRUE)
comm.print(tmp)
tmp <- comm.all(a, na.rm = TRUE)
comm.print(tmp)

### Finish.
finalize()
