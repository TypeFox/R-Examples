### This is a famous example from snow package.
### One can source this file into R interactive model or run it by the command
### SHELL> Rscript --vanilla 02_spmd.M_spmd.r

library(parallel)
cl <- makeCluster(2)

splitRows <- function (x, ncl){
  lapply(splitIndices(nrow(x), ncl), function(i) x[i, , drop = FALSE])
}
parMM <- function (cl, A, B){
  do.call(rbind, clusterApply(cl, splitRows(A, length(cl)), get("%*%"), B))
}

set.seed(123)
A <- matrix(rnorm(1000000), 1000)
system.time(replicate(10, A %*% A))
system.time(replicate(10, parMM(cl, A, A)))

stopCluster(cl)
