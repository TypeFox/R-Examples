# Initialize MPI
suppressPackageStartupMessages(library(pbdDEMO, quietly = TRUE))
init.grid()

# Number of rows/columns
n <- 250
p <- 50

# Generate a global matrix on processor 0, then distribute
comm.set.seed(1234, diff = TRUE)

# ICTXT=0
dx <- diag(1, 10, 10, type="ddmatrix", bldim=2)

comm.cat("\nthe diagonal matrix on context 0\n", quiet=TRUE)
print(dx)
comm.print(submatrix(dx), all.rank=TRUE)

# ICTXT=2
dx <- redistribute(dx, ICTXT=2)

comm.cat("\n\n\nthe diagonal matrix on context 2\n", quiet=TRUE)
print(dx)
comm.print(submatrix(dx), all.rank=TRUE)

# undistributed
comm.cat("\n\n\nthe matrix as a global, non-distributed matrix\n", quiet=TRUE)
x <- as.matrix(dx)

comm.print(x, quiet=TRUE)

finalize()
