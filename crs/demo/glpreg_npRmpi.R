## Demo file for running parallel npglpreg. Setting mpi=TRUE in
## npglpreg will call npRmpi rather than np for computing the
## estimator. Other than that running this example is identical to
## running npRmpi code. Note you _must_ have .Rprofile in your current
## directory (rename Rprofile from the npRmpi inst directory to
## .Rprofile then follow instructions for running parallel R jobs on
## your system, e.g. openmpirun -n 4 R CMD BATCH glpreg_npRmpi.R)

rm(list=ls())

mpi.bcast.cmd(np.mpi.initialize(),
              caller.execute=TRUE)

set.seed(42)
n <-  250

degree.max <- 20

mpi.bcast.cmd(library(crs),
              caller.execute=TRUE)

x1 <- runif(n)
x2 <- runif(n)
dgp <- cos(8*pi*x1)
y <- dgp+rnorm(n,sd=0.1)
  
mpi.bcast.Robj2slave(y)
mpi.bcast.Robj2slave(x1)
mpi.bcast.Robj2slave(x2)
mpi.bcast.Robj2slave(degree.max)

mpi.bcast.cmd(model.glp <- npglpreg(y~x1+x2,degree.max=degree.max,mpi=TRUE),
              caller.execute=TRUE)

summary(model.glp)

mpi.bcast.cmd(mpi.quit(),
              caller.execute=TRUE)
