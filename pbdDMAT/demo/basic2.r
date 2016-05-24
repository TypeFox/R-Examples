### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

# Initialize process grid
library(pbdDMAT, quiet=T)

if(comm.size() != 2)
  comm.stop("Exactly 2 processors are required for this demo.")

init.grid()

# First we generate the matrix on process 0 and then distribute it to
  # all the others
{
  if (comm.rank()==0)
    A <- 1:100
  else
    A <- NULL
}

dA <- as.ddmatrix(A, bldim=2) # distribute

# Notice how the vector is stored on a 1x2 grid
print(dA)
comm.print(head(submatrix(dA)))

# Let's convert it to a 2x1 grid; context 2 is your default (nprocs)x1 grid
# Notice that the data is distributed differently
newA <- redistribute(dA, bldim=bldim(dA), ICTXT=2)
print(newA)
comm.print(head(submatrix(newA)))

# Be careful when redistributing things across different contexts
# When performing operations involving two distributed matrices, such
  # as A %*% B, they must be distributed on the same context and
  # with the same blocking factor!

# Finish
finalize()
