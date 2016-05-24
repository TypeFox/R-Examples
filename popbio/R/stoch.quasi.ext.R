stoch.quasi.ext<-function(matrices, n0, Nx, tmax=50, maxruns=10, nreps=5000, prob=NULL, sumweight=NULL, verbose=TRUE)
{
   if(is.list(matrices)){matrices<-matrix(unlist(matrices), ncol=length(matrices))}
   x<-length(n0)
   if(is.null(sumweight)){sumweight=rep(1,x)}
   y<-dim(matrices)[2]
   ext<-matrix(numeric(maxruns*tmax), ncol=maxruns)   
   for(h in 1:maxruns)
   {
      if(verbose)
      {
         print(paste("Calculating extinction probability for run", h), quote=FALSE)
      }
      prob.ext<-numeric(tmax)
      for(i in 1:nreps)
      {
         n<-n0
         for( t in 1:tmax)
         {
            col<-sample(1:y, 1, prob= prob)
            A<-matrix(matrices[,col], nrow=x )
            n<-A %*% n
            N<-sum(sumweight*round(n))
            if(N<Nx)
            {
               prob.ext[t]<-prob.ext[t]+1
               break
            }
         }      
      }
     prob.ext<-cumsum(prob.ext/nreps)
     ext[,h]<-prob.ext   
   }
ext
}

