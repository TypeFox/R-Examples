### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

# Initialize process grid
library(pbdDMAT, quiet=T)

if(comm.size() != 2)
  comm.stop("Exactly 2 processors are required for this demo.")

init.grid()

# First we generate the matrix on process 0 and then distribute it to
  # all the others
if (comm.rank()==0){
  A <- 1:100
  nas <- sample(A, replace=F, size=5)
  A[which(A %in% nas)] <- NA
  A <- matrix(A, 10, 10)
  print(A)
} else {
  A <- NULL
}

dA <- as.ddmatrix(A, bldim=2) # distribute

comm.print(dA, all.rank=T)

print(dA)
comm.print(head(submatrix(dA)))

# Throw away rows with NA's in them
dA <- na.exclude(dA)

comm.print(attr(dA, "na.action"))

print(dA)

# Check against R to make sure we did it right
ourA <- as.matrix(dA)#, proc.dest=0)

if (comm.rank()==0){
  trueA <- na.exclude(A)
  
  print(trueA)
  print(ourA)
  
  cat(sprintf("\nDid we exclude NA's correctly?\n"))
  all(trueA==ourA) # we could use all.equal(), but there's no 
    # guarantee that the 'na.action' attributes will agree on order
    # and so we just compare the elements of the resulting matrices.
}

# Finish
finalize()
