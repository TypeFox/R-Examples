##create a matrix
mat <- cbind(c(1,0,0),c(0,2,1),c(0,1,2))
##define the number of blocks and block sizes
block.sizes <- c(1,2)
n.blocks <- length(block.sizes)

##Compute the Cholesky factor
R <- makeCholBlock(mat, n.blocks, block.sizes)
##and the matrix inverse
i.mat <- invCholBlock(R, n.blocks, block.sizes)
##compare to the alternative
i.mat-solve(mat)
\dontshow{
  if( max(abs(R-chol(mat))) > 1e-10 ){
    stop("makeCholBlock: Results not equal")
  }
  if( max(abs(i.mat-solve(mat))) > 1e-10 ){
    stop("invCholBlock: Results not equal")
  }
}
##define a B vector
B <- c(1,2,3)
##solve the equation system (we need n.x since B is not a matrix)
x1 <- solveTriBlock(R, B, n.blocks, block.sizes, tr=TRUE)
x2 <- solveTriBlock(R, x1, n.blocks, block.sizes, tr=FALSE)
print(x2)
##compare to the alternative
print(solve(mat,B))
range(x2-solve(mat,B))
\dontshow{
  if( max(abs(x2-solve(mat,B))) > 1e-10 ){
    stop("solveTriBlock: Results not equal")
  }
}
##compute the quadratic form B'*i.mat*B
norm2(x1)
##compare to the alternative
t(B) %*% i.mat %*% B
\dontshow{
  if( abs(norm2(x1) - t(B) %*% i.mat %*% B) > 1e-10 ){
    stop("solveTriBlock: Results not equal for norm")
  }
}

