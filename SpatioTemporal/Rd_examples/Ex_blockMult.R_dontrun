#create a matrix
mat <- cbind(c(1,0,0), c(0,2,1), c(0,1,2))
#define the number of blocks and block sizes
block.sizes <- c(1,2)
n.blocks <- length(block.sizes)

#define a X vector
X <- matrix(c(1,2,3,1,1,1), 3, 2)

#compute mat %*% X
blockMult(mat, X, n.blocks, block.sizes)
#or the old fashioned way
mat %*% X
\dontshow{
  if( max(abs(blockMult(mat, X, n.blocks, block.sizes) -
              (mat %*% X))) > 1e-13 ){
    stop("blockMult: Results not equal")
  }
}

