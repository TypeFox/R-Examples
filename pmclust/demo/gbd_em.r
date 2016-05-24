### Setup mpi environment.
suppressMessages(library(pmclust, quietly = TRUE))
comm.set.seed(123, diff = TRUE)

### Generate an example data.
N.allgbds <- rep(5000, comm.size())
N.gbd <- 5000
N.K.gbd <- c(2000, 3000)
N <- 5000 * comm.size()
p <- 2
K <- 2
data.gbd <- generate.basic(N.allgbds, N.gbd, N.K.gbd, N, p, K)

### Run clustering.
PARAM.org <- set.global.gbd(K = K, X.gbd = data.gbd$X.spmd)
PARAM.org <- initial.em(PARAM.org)
PARAM.new <- em.step(PARAM.org)
em.update.class()
mb.print(PARAM.new, .pmclustEnv$CHECK)

### Get results.
N.CLASS <- get.N.CLASS(K)
comm.cat("# of class:", N.CLASS, "\n")

### Print run time and quit Rmpi.
finalize()
