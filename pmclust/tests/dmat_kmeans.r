### Setup mpi environment.
suppressMessages(library(pbdDMAT, quietly = TRUE))
suppressMessages(library(pmclust, quietly = TRUE))
init.grid()
comm.set.seed(123, diff = TRUE)

### Generate an example data.
N.allspmds <- rep(5000, comm.size())
N.spmd <- 5000
N.K.spmd <- c(2000, 3000)
N <- 5000 * comm.size()
p <- 2
K <- 2
data.spmd <- generate.basic(N.allspmds, N.spmd, N.K.spmd, N, p, K)
X.dmat <- as.dmat(data.spmd$X.spmd)

### Run clustering.
PARAM.org <- set.global.dmat(K = K)
PARAM.org <- initial.center.dmat(PARAM.org)
PARAM.new <- kmeans.step.dmat(PARAM.org)
kmeans.update.class.dmat()
mb.print(PARAM.new, .pmclustEnv$CHECK)

### Get results.
N.CLASS <- get.N.CLASS.dmat(K)
comm.cat("# of class:", N.CLASS, "\n")

### Print run time and quit Rmpi.
finalize()
