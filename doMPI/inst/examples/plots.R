# This example shows how to write a parallel program where the results
# are a sequence of plots, stored as PNG files.  It's the sort of
# example that you could run as a sequence batch jobs, but I think this
# is much simpler, as well as more portable.
#
# Note: This example can fail on Mac OS X.  I get the following
# error message when starting remote workers with orterun:
#
#  "On-demand launch of the Window Server is allowed for root user only."
#
# Presumably this is due to the use of the PNG device, but I haven't
# tracked it down.

suppressMessages(library(doMPI))

# Create and register an MPI cluster
cl <- startMPIcluster()
registerDoMPI(cl)

# Initialize variables
trials <- 10
n <- nrow(iris)
leaveout <- 2

# Define chunkSize so that each cluster worker gets a single "task chunk"
chunkSize <- ceiling(trials / getDoParWorkers())
mpiopts <- list(chunkSize=chunkSize)

# Define a .combine function that throws away the "results"
trash <- function(...) NULL

# Create the PNG files in parallel
foreach(i=icount(trials), .combine=trash, .multicombine=TRUE,
        .packages='randomForest', .options.mpi=mpiopts) %dopar% {
  d <- iris[sample(n, n - leaveout),]
  rf <- randomForest(Species~., data=d, proximity=TRUE)
  png(filename=sprintf('MDSplot_%d.png', i))
  MDSplot(rf, d$Species)
  dev.off()
  NULL
}

# Shutdown the cluster and quit
closeCluster(cl)
mpi.quit()
