### Setup environment.
suppressMessages(library(pmclust, quietly = TRUE))
comm.set.seed(123, diff = TRUE)

### Generate an example data.
N <- 5000
p <- 2
K <- 2
data.spmd <- generate.MixSim(N, p, K, MaxOmega = 0.001)
X.spmd <- data.spmd$X.spmd

### Run clustering.
PARAM.org <- set.global(K = K)          # Set global storages.
 PARAM.org <- initial.em(PARAM.org)    # One initial.
# PARAM.org <- initial.RndEM(PARAM.org)   # Ten initials by default.
PARAM.new <- apecma.step(PARAM.org)     # Run APECMa.
em.update.class()                       # Get classification.

### Get results.
N.CLASS <- get.N.CLASS(K)
comm.cat("# of class:", N.CLASS, "\n", quiet = TRUE)
comm.cat("# of class (true):", data.spmd$N.CLASS, "\n", quiet = TRUE)

### Quit.
finalize()

