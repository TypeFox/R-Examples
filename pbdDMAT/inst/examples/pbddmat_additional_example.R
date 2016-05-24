library(pbdDMAT, quiet=T)

# Demonstrating generation of matrices without communication
# (what is needed is generated locally), calculating the 
# singular value decomposition in parallel, and timing parallel
# operations. Several of these functions are arguably useful but 
# not included in any existing pbd package, though they are slated
# for a future package.

# This script generates a random normal matrix of specified size
# and with specified parameters independently on the processors 
# (so that there is no data communication in matrix generation)
# using independent seeds via the rlecuyer package. Then a 
# principal components analysis is performed via prcomp() which 
# uses svd() --- powered by ScaLAPACK.

# This analysis is timed and the minimum, maximum, and mean run 
# times across all processors is reported.

# New users are encouraged to experiment with different processor
# grid shapes (vid nprow= and npcol= in init.grid()), blocking 
# factors, and matrix sizes.

# The singular value decomposition is an extremely computationally 
# intensive operation, being O(n^3) in serial for a square matrix 
# of dimension nxn. An analysis of SVD complexity can be found in 
# Lecture 31, pages 239-244 of "Numerical Linear Algebra" by 
# Lloyd N. Trefethen and David Bau, III

# ---------------------------------------------
# Setup
# ---------------------------------------------

init.grid() 

# Use rlecuyer to generate independent random seeds on the processors
g.set.seed(diff=T)

# Number of rows and columns to generate
nrows <- 1e3
ncols <- 1e3

mn <- 10
sdd <- 100

# ScaLAPACK blocking dimension
.BLDIM <- c(8, 8)

# ---------------------------------------------
# Analysis
# ---------------------------------------------

# timing function
timer <- function(timed)
{
  ltime <- system.time(timed)[3]
  mintime <- allreduce(ltime, op='min')
  maxtime <- allreduce(ltime, op='max')
  
  meantime <- allreduce(ltime, op='sum') / comm.size()
  
  return( c(min=mintime, mean=meantime, max=maxtime) )
}

# generated
time_data <- timer({
  x <- ddmatrix("rnorm", nrow=nrows, ncol=ncols, mean=mn, sd=sdd)
})

barrier()

time_pca <- timer({
  pca_x <- prcomp(x, scale=TRUE)
})

barrier()

comm.print( 
  rbind(
    DataGeneration=time_data, 
    PCA=time_pca, 
    Total=time_pca+time_data), 
  quiet=T)

finalize()
