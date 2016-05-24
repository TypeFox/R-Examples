# Rdsm code to find max burst in a time series; finds the block in the
# input vector having maximal average value, among blocks of size k

# for example, in (5,7,6,20,4,14,11,12,15,17), the maximal block,
# (15,17), starts at element 9 and has the average value 16

# requires the zoo package, which performs various simply time series
# operations

# arguments:

#    x:  data vector
#    k:  block size
#    mas:  scratch space, shared, 1 x (length(x)-1)
#    rslts:  2-tuple showing the maximum burst value, and 
#            where it starts; shared, 1 x 2

maxburst <- function(x,k,mas,rslts) {
   require(Rdsm)
   require(zoo)
   # determine this thread's chunk of x
   n <- length(x)
   myidxs <- getidxs(n-k+1)
   myfirst <- myidxs[1]
   mylast <- myidxs[length(myidxs)]
   mas[1,myfirst:mylast] <- rollmean(x[myfirst:(mylast+k-1)],k)
   barr()  # make sure all threads have written to mas
   # one thread does wrapup, might as well be thread 1
   if (myinfo$id == 1) {
      rslts[1,1] <- which.max(mas[,])
      rslts[1,2] <- mas[1,rslts[1,1]]
   }
}

test <- function(cls) {
   require(Rdsm)
   mgrinit(cls)
   mgrmakevar(cls,"mas",1,9)
   mgrmakevar(cls,"rslts",1,2)
   x <<- c(5,7,6,20,4,14,11,12,15,17)
   clusterExport(cls,"maxburst")
   clusterExport(cls,"x")
   clusterEvalQ(cls,maxburst(x,2,mas,rslts))
   print(rslts[,])  # not print(rslts)!
}
