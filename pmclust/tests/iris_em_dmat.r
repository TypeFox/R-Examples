# File name: iris_dmat.r
# Run: mpiexec -np 4 Rscript iris_dmat.r

rm(list = ls())                                       # Clean environment
suppressMessages(library(pbdDMAT, quietly = TRUE))                      # Load library
init()
init.grid()

### Load data
X <- as.matrix(iris[, -5])                            # Dimension 150 by 4
X.cid <- as.numeric(iris[, 5])                        # True id

### Convert to ddmatrix
X.dmat <- as.ddmatrix(X)

### Standardized
X.std <- scale(X.dmat)

### Clustering
suppressMessages(library(pmclust, quietly = TRUE))
comm.set.seed(123, diff = TRUE)

X.dmat <- X.std
PARAM.org <- set.global.dmat(K = 3)                   # Preset storage
.pmclustEnv$CONTROL$debug <- 1                        # Disable debug messages
PARAM.org <- initial.em.dmat(PARAM.org)
# PARAM.org <- initial.RndEM.dmat(PARAM.org)
PARAM.mbc <- em.step.dmat(PARAM.org)                  # model-based
em.update.class.dmat()
mb.print(PARAM.mbc, .pmclustEnv$CHECK)

### Get results.
N.CLASS <- get.N.CLASS.dmat(K = 3)
comm.cat("# of class:", N.CLASS, "\n")

### Finish
finalize()
