library(pbdBASE, quiet=TRUE)

init.grid()

nrows <- ncols <- 10
BL <- 2

# Generate data on process 0, then distribute to the others
{
if (comm.rank()==0)
   x <- matrix(rnorm(n=nrows*ncols, mean=10, sd=100), nrow=nrows, ncol=ncols)
else
  x <- NULL
}
  
dx <- as.ddmatrix(x=x, bldim=BL)

print(dx)
comm.cat("\n", quiet=T) # put a blank line for easier reading of the output

comm.print(submatrix(dx))
comm.cat("\n", quiet=T)

comm.print(dx)
 
dx[1,1] <- NA # insertion indices are global
comm.print(submatrix(dx)[1,1], all.rank=T) # see?
comm.cat("\n", quiet=T)

comm.print(dim(dx))
dx <- dx[, -2]
comm.print(dim(dx))
comm.cat("\n", quiet=T)

nona <- na.exclude(dx)

# convert back
nona <- as.matrix(nona, proc.dest=0)

# compare our results with R --- notice the syntax is essentially identical
if (comm.rank()==0){
  x[1,1] <- NA
  x <- x[, -2]
  r_nona <- na.exclude(x)
  
  all.equal(r_nona, nona)
}

finalize()
