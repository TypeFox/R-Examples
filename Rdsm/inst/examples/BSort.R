# bucket sort with sampling

# vector x is broken into chunks according to cut points; this implies
# that all the numbers in chunks i are <= all those in chunk j > i; thus
# each chunk can be sorted and then placed into its proper place in x

# the cuts are obtained by first sorting a sample of x and then
# computing r - 1 quantiles, where r is the number of threads

# arguments:
# 
#    x:  vector to be sorted; shared; sorted in place
#    counts:  intermediate result; shared, length = length(cls)
#    samp:  intermediate result; shared; length = nsamp
#    nsamp:  number of elements of x to sample

bsort <- function(x,counts,samp,nsamp=1000) {
   me <- myinfo$id
   # make local copy of x to avoid cache coherency overhead
   tmpx <- x[1,]
   if (me == 1) {  # sample to get quantiles
      samp[1,] <- sort(tmpx[sample(1:length(tmpx),nsamp,replace=F)])
   }
   barr()
   # determine my interval
   r <- myinfo$nwrkrs
   k <- floor(nsamp / r)
   if (me > 1) mylo <- samp[1,(me-1) * k]
   if (me < r) myhi <- samp[1,me * k]
   # get my chunk and sort it
   if (me == 1) myx <- tmpx[tmpx <= myhi] else 
      if (me == r) myx <- tmpx[tmpx > mylo] else
         myx <- tmpx[ tmpx > mylo & tmpx <= myhi]
   myx <- sort(myx)
   # need to decide where in x to place myx
   lx <- length(myx)
   counts[1,me] <- lx
   barr()
   if (me == 1) counts[1,] <- cumsum(counts[1,])
   barr()
   # place my sorted chunk back in x
   if (me == 1) x[1,1:lx] <- myx else {
      start <- counts[1,me-1] + 1
      fin <- start + lx - 1
      x[1,start:fin] <- myx
   }
}

test <- function(cls,barrback=F) {
   require(parallel)
   mgrinit(cls,barrback=barrback)
   mgrmakevar(cls,"a",1,25,fs=barrback)
   mgrmakevar(cls,"counts",1,length(cls),fs=barrback)
   mgrmakevar(cls,"smp",1,10,fs=barrback)
   a[1,] <- runif(25)
   print(a[1,])
   clusterExport(cls,"bsort")
   clusterEvalQ(cls,bsort(a,counts,smp,nsamp=10))
   print(a[1,])
}

