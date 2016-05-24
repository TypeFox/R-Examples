# This is a simple parallel random forest example, using Andy Liaw's
# randomForest package.

suppressMessages(library(doMPI))
suppressMessages(library(randomForest))

# Create and register an MPI cluster
cl <- startMPIcluster()
registerDoMPI(cl)

# Define a parallel randomForest function
rforest <- function(x, y=NULL, xtest=NULL, ytest=NULL, ntree=500, ...) {
  initWorkers <- function() library(randomForest)
  opts <- list(initEnvir=initWorkers)

  foreach(i=idiv(ntree, chunks=getDoParWorkers()),
          .combine='combine', .multicombine=TRUE, .inorder=FALSE,
          .options.mpi=opts) %dopar% {
    randomForest:::randomForest.default(x, y, xtest, ytest, ntree=i, ...)
  }
}

# Create a matrix and factor as input
m <- 200; n <- 120
x <- matrix(rnorm(m * n), m, n)
y <- gl(10, m/10)

# Execute rforest and then print the resulting model object
rfit <- rforest(x, y)
print(rfit)

# Shutdown the cluster and quit
closeCluster(cl)
mpi.quit()
