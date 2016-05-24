library(doRNG)
library(doParallel)
registerDoParallel(cores=2)

# single %dorng% loops are reproducible
r1 <- foreach(i=1:4, .options.RNG=1234) %dorng% { runif(1) }
r2 <- foreach(i=1:4, .options.RNG=1234) %dorng% { runif(1) }
identical(r1, r2)
# the sequence os RNG seed is stored as an attribute
attr(r1, 'rng')

# sequences of %dorng% loops are reproducible
set.seed(1234)
s1 <- foreach(i=1:4) %dorng% { runif(1) }
s2 <- foreach(i=1:4) %dorng% { runif(1) }
# two consecutive (unseed) %dorng% loops are not identical
identical(s1, s2)

# But the whole sequence of loops is reproducible
set.seed(1234)
s1.2 <- foreach(i=1:4) %dorng% { runif(1) }
s2.2 <- foreach(i=1:4) %dorng% { runif(1) }
identical(s1, s1.2) && identical(s2, s2.2)
# it gives the same result as with .options.RNG
identical(r1, s1)

# Works with SNOW-like and MPI clusters
# SNOW-like cluster
cl <- makeCluster(2)
registerDoParallel(cl)

s1 <- foreach(i=1:4, .options.RNG=1234) %dorng% { runif(1) }
s2 <- foreach(i=1:4, .options.RNG=1234) %dorng% { runif(1) }
identical(s1, s2)

stopCluster(cl)
registerDoSEQ()

# MPI cluster
library(doMPI)
cl <- startMPIcluster(2)
registerDoMPI(cl)

s1 <- foreach(i=1:4, .options.RNG=1234) %dorng% { runif(1) }
s2 <- foreach(i=1:4, .options.RNG=1234) %dorng% { runif(1) }
identical(s1, s2)

closeCluster(cl)
registerDoSEQ()
