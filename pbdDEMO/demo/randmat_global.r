### SHELL> mpiexec -np 4 Rscript --vanilla [...].r

# Initialize MPI
library(pbdDEMO, quietly = TRUE)
init.grid()

# Number of rows/columns
n <- 250
p <- 50

# Generate a global matrix on processor 0, then distribute
comm.set.seed(1234, diff = TRUE)
if (comm.rank()==0){
  x1 <- matrix(rnorm(n*p), nrow=n, ncol=p)
} else {
  x1 <- NULL
}

dx1 <- as.ddmatrix(x1)

# Generate _the same_ global matrix on all processors, then distribute
# not useful in application, but handy for testing
comm.set.seed(1234, diff = FALSE)
x2 <- matrix(rnorm(n*p), nrow=n, ncol=p)

dx2 <- as.ddmatrix(x2)


# See that these two method are equivalent
comm.cat("Are these the same distributed matrix?\n", quiet=T)
test <- all.equal(dx1, dx2)
if (is.logical(test)){
  comm.cat("YES!\n\n", quiet=T)
} else {
  comm.cat("No...\n", quiet=T)
  comm.print(test, quiet=T)
}



# Generate locally only what is needed.
# This will produce a different matrix because of the block cyclic 
# distribution.
comm.set.seed(1234, diff = TRUE)
dx3 <- ddmatrix("rnorm", nrow=n, ncol=p, bldim=4)

print(dx3)

finalize()
