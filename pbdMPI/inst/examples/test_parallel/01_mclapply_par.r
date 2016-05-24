### This is a famous example from multicore package.
### One can source this file into R interactive model or run it by the command
### SHELL> Rscript --vanilla 01_mclapply_spmd.r

library(parallel)

system.time(
  ret <- unlist(mclapply(1:32, function(x) sum(rnorm(1e7))))
)
