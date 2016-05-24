
# matrix multiplication; the product u %*% v is computed on cls, and
# stored in w; u, v and w are big.matrix objects

# the result will be in w, rather than being returned by the call

# simpler test:  direct launch of threads, no wrapper
test <- function(cls) {
   require(parallel)
   require(Rdsm)
   mgrinit(cls)  # start Rdsm
   # create the shared variables, and put data in them
   mgrmakevar(cls,"a",6,2)
   mgrmakevar(cls,"b",2,6)
   mgrmakevar(cls,"c",6,6)
   a[,] <- 1:12
   b[,] <- 1  
   # define the function to be executed by each thread
   mmul1thread <- function(u,v,w) {
      require("parallel")
      # determine which rows this thread will process
      myidxs <- splitIndices(nrow(u),myinfo$nwrkrs)[[myinfo$id]]
      # now process those rows
      w[myidxs,] <- u[myidxs,] %*% v[,]
   }
   clusterExport(cls,"mmul1thread",envir=environment())
   # have the threads do their work
   clusterEvalQ(cls,mmul1thread(a,b,c))
   print(c[,])  # check result
}

# more general test, differing from test(), with  a wrapper function
# mmul() (see below); also, this test uses filebacked variables
test1 <- function(cls) {
   require(parallel)
   require(Rdsm)
   mgrinit(cls)
   mgrmakevar(cls,"a",6,2,fs=TRUE)
   mgrmakevar(cls,"b",2,6,fs=TRUE)
   mgrmakevar(cls,"c",6,6,fs=TRUE)
   a[,] <- 1:12
   b[,] <- 1  # all 1s
   mmul("a","b","c",cls)
   print(c[,])
}

# matrices specified by quoted name
mmul <- function(uname,vname,wname,cls) {
   require(parallel)
   clusterCall(cls,mmul1thread,uname,vname,wname)
}

# code to be run by a thread; u, v, w can be specified in any 
# form--bigmatrix object, bigmatrix decriptor, quoted bigmatrix 
# name--resolved by getmatrix()
mmul1thread<- function(u,v,w) {
   require(parallel)
   require(Rdsm)
   u <- getmatrix(u)
   v <- getmatrix(v)
   w <- getmatrix(w)
   myidxs <- splitIndices(nrow(u),myinfo$nwrkrs)[[myinfo$id]]
   w[myidxs,] <- u[myidxs,] %*% v[,]
   invisible(0)  # don't do expensive return of result
}

