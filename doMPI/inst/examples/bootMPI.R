# This example comes from an old version of the Wikipedia article on
# bootstrapping, which actually included R code.  The task being executed
# in the foreach loop executes very quickly, but it does a lot of trials.
# That makes it a good example of the need for the chunkSize option.

suppressMessages(library(doMPI))

# Create and register an MPI cluster
cl <- startMPIcluster()
registerDoMPI(cl)

# Initialize variables
x <- iris[which(iris[,5] != "setosa"), c(1,5)]
trials <- 10000

# Define chunkSize so that each cluster worker gets a single "task chunk"
chunkSize <- ceiling(trials / getDoParWorkers())
mpiopts <- list(chunkSize=chunkSize)

# Perform the bootstrapping in parallel
r <- foreach(icount(trials), .combine='cbind', .options.mpi=mpiopts) %dopar% {
  ind <- sample(100, 100, replace=TRUE)
  result1 <- glm(x[ind,2]~x[ind,1], family=binomial(logit))
  structure(coefficients(result1), names=NULL)
}

# Print the resulting matrix
print(r)

# Shutdown the cluster and quit
closeCluster(cl)
mpi.quit()
