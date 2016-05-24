# k-means clustering on the data matrix x, with k clusters and ni
# iterations; final cluster centroids placed in cntrds

# initial centroids taken to be k randomly chosen rows of x; if a
# cluster becomes empty, its new centroid will be a random row of
# x

# arguments:
#    x:  data matrix x; shared
#    k:  number of clusters
#    ni:  number of iterations
#    cntrds:  centroids matrix; row i is centroid i; shared, k by ncol(x)
#    cinit:  optional initial values for the centroids; k by ncol(x)
#    sums:  scratch matrix; sums[j,] contains the count 
#       and sum for cluster j; shared, k by 1+ncol(x)
#    lck:  lock variable; shared

# the result will be in cntrds

kmeans <- function(x,k,ni,cntrds,sums,lck,cinit=NULL) {
   require(parallel)
   require(pdist)
   nx <- nrow(x)
   # get my assigned portion of x
   # myidxs <- splitIndices(nx,myinfo$nwrkrs)[[myinfo$id]]
   myidxs <- getidxs(nx)
   myx <- x[myidxs,]  
   # random initial centroids if none specified
   if (is.null(cinit)) {
      if (myinfo$id == 1) 
         cntrds[,] <- x[sample(1:nx,k,replace=F),]
      barr()
   } else cntrds[,] <- cinit

   # mysum()sums the rows in myx corresponding to the indices idxs; we
   # also produce a count of those rows
   mysum <- function(idxs,myx) {
      c(length(idxs),colSums(myx[idxs,,drop=F]))
   }
   for (i in 1:ni) {  # ni iterations
      # node 1 is sometimes asked to do some "housekeeping"
      if (myinfo$id == 1) {
         sums[] <- 0
      }
      barr()  # other nodes wait for node 1 to do its work
      # find distances from my rows of x to the centroids, then
      # find which centroid is closest to each such row
      dsts <- matrix(pdist(myx,cntrds[,])@dist,ncol=nrow(myx))
      nrst <- apply(dsts,2,which.min)
      # nrst[i] contains the index of the nearest centroid to row i in
      # myx
      tmp <- tapply(1:nrow(myx),nrst,mysum,myx)
      # in the above, we gather the observations in myx whose closest
      # centroid is centroid j, and find their sum, placing it in
      # tmp[j]; the latter will also have the count of such observations
      # in its leading component 
      # next, we need to add that to sums[j,], as an atomic operation
      rdsmlock(lck)
      # the j values in tmp will be strings, so convert
      for (jn in names(tmp)) {
         j <- as.integer(jn)
         sums[j,] <- sums[j,] + tmp[[jn]]
      }
      rdsmunlock(lck)
      barr()  # wait from sums[,] to be ready
      if (myinfo$id == 1) {
         # update centroids, using a random data point if a cluster
         # becomes empty
         for (j in 1:k) {
           # update centroid for cluster j
           if (sums[j,1] > 0) {
              cntrds[j,] <- sums[j,-1] / sums[j,1] 
            } else cntrds[j] <<- x[sample(1:nx,1),]
         }
      }
   }
   0  # don't do expensive return of result
}

test <- function(cls,boost=F) {
   require(parallel)
   mgrinit(cls,boost)
   mgrmakevar(cls,"x",12,2)
   mgrmakevar(cls,"cntrds",2,2)
   mgrmakevar(cls,"sms",2,3)
   mgrmakelock(cls,"lck",boost)
   x[,] <- matrix(sample(1:100,24,replace=TRUE),ncol=2)
   clusterExport(cls,"kmeans")
   # clusterEvalQ(cls,debug(kmeans))
   if (boost) {
      clusterEvalQ(cls,kmeans(x,2,1,cntrds,sms,lck,
         cinit=rbind(c(5,5),c(15,15))))
   } else {
      clusterEvalQ(cls,kmeans(x,2,1,cntrds,sms,"lck",
         cinit=rbind(c(5,5),c(15,15))))
   }
   print(cntrds[,])
}

test1 <- function(cls) {
   require(parallel)
   mgrinit(cls,boost=TRUE,barrback=TRUE)
   mgrmakevar(cls,"x",10000,3)
   mgrmakevar(cls,"cntrds",3,3)
   mgrmakevar(cls,"sms",3,4)
   mgrmakelock(cls,"lck",boost=TRUE)
   x[,] <- matrix(rnorm(30000),ncol=3)
   ri <- sample(1:10000,3000)
   x[ri,1] <- x[ri,1] + 5
   ri <- sample(1:10000,3000)
   x[ri,2] <- x[ri,2] + 5
   clusterExport(cls,"kmeans")
   clusterEvalQ(cls,kmeans(x,3,50,cntrds,sms,lck))
   print(cntrds[,])
}
