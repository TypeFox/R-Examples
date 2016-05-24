### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

# Initialize process grid
library(pbdDMAT, quiet=T)

if(comm.size() != 2)
  comm.stop("Exactly 2 processors are required for this demo.")

init.grid()

# Setup for the remainder
set.seed(25)
M <- N <- 6 # number of rows/cols
BL <- 2 # blocking factor --- passing single value BL assumes BLxBL blocking

# First we generate the matrix on process 0 and then distribute it to
  # all the others
{
  if (comm.rank()==0)
    A <- matrix(rnorm(M * N, mean = 100, sd = 10), nrow = M, ncol = N)
  else
    A <- NULL
}

dA <- as.ddmatrix(A, bldim=BL) # distribute

# Print some information about the matrix
print(dA)

# Set the first column to zero of the global matrix
dA[, 1] <- 0

print(dA)

newA <- as.matrix(dA)
comm.print(newA)

# Get a global vector of the first column, stored only on process 0
# Whenever a non-distributed matrix is owned by only one processor, 
  # the other processors store 'NA' for that object
first.col <- dA[,1] # subsetting as a distributed object
fc.gl <- as.vector(first.col, proc.dest=0) # convert to global on process 0
comm.print(fc.gl, all.rank=T)


# Finish
finalize()
