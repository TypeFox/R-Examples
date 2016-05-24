# This example shows how to use the initEnvir option to initialize
# the cluster workers.  This is necessary when the workers need to
# use an object that can't be simply exported, such as connection
# objects.
#
# The function specified by the initEnvir option is passed the
# worker's execution environment, plus any arguments specified
# via the initArgs option.
#
# In this example, I demonstrate how initEnvir and finalEnvir
# can be used to create a file for each of the workers that
# can be accessed via the variable "out".  This can be used
# to write task results from the workers without opening and
# closing the file for each task, and without conflicts
# between the different workers.  This technique is particularly
# useful if the job is executed from an NFS directory, so that
# all of the files are accessible by the user at the end of
# the job.

suppressMessages(library(doMPI))

# Create and register an MPI cluster
cl <- startMPIcluster()
registerDoMPI(cl)

# Create an output file that the workers can write to using
# the variable "out"
initEnvir <- function(envir, template, comm) {
  rank <- mpi.comm.rank(comm)
  fname <- sprintf(template, rank)
  envir$out <- file(fname, 'w')
  cat(sprintf('This is the output file for worker %d\n', rank), file=envir$out)
}

# Pass the template of the output file name to the initEnvir function
initArgs <- list(template='worker_%d.out', comm=cl$comm)

# Close the output file at the end of the foreach loop
finalEnvir <- function(envir) {
  close(envir$out)
}

# Define a .combine function that throws away it's arguments
trash <- function(...) NULL

# Define the doMPI-specific options
mpiopts <- list(initEnvir=initEnvir, initArgs=initArgs, finalEnvir=finalEnvir)

# Specify the doMPI-specific opitons with ".options.mpi"
# Note that "trash" can be called with multiple arguments
foreach(i=1:10, .options.mpi=mpiopts,
        .combine=trash, .multicombine=TRUE) %dopar% {
  cat(sprintf('This is task number %d\n', i), file=out)
}

# Shutdown the cluster and quit
closeCluster(cl)
mpi.quit()
