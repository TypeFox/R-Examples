
# finds distances between all possible pairs of rows in the matrix x and
# rows in the matrix y, as with pdist() but in parallel

# arguments:
#    x,y: input matrices
#    dout:  output matrix, shared

rdsmpdist <- function(x,y,dout) {
   require(pdist)
   nx <- nrow(x)
   myidxs <-  getidxs(nx)
   dout[myidxs,] <- as.matrix(pdist(x[myidxs,],y[,]))
   invisible(0)  # don't do expensive return of result
}

test <- function(cls) {
   require(parallel)
   require(pdist)
   mgrinit(cls)
   mgrmakevar(cls,"a",4,2)
   mgrmakevar(cls,"b",2,4)
   mgrmakevar(cls,"c",4,2)
   a[,] <- 1:8
   b[,] <- 1  # all 1s
   clusterExport(cls,"rdsmpdist")
   clusterEvalQ(cls,rdsmpdist(a,b,c))
   print(c[,])
   print(as.matrix(pdist(a[,],b[,])))
}

