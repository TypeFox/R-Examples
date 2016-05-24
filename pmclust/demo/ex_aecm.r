### Setup mpi environment.
suppressMessages(library(pmclust, quietly = TRUE))
comm.set.seed(123, diff = TRUE)

### Generate an example data.
N.allspmds <- rep(5000, comm.size())
N.spmd <- 5000
N.K.spmd <- c(2000, 3000)
N <- 5000 * comm.size()
p <- 2
K <- 2
data.spmd <- generate.basic(N.allspmds, N.spmd, N.K.spmd, N, p, K)
X.spmd <- data.spmd$X.spmd

### Run clustering.
PARAM.org <- set.global(K = K)
PARAM.org <- initial.em(PARAM.org)
PARAM.new <- aecm.step(PARAM.org)
em.update.class()
mb.print(PARAM.new, .pmclustEnv$CHECK)

### Get results.
N.CLASS <- get.N.CLASS(K)
comm.cat("# of class:", N.CLASS, "\n")

### Print run time and quit Rmpi.
finalize()
