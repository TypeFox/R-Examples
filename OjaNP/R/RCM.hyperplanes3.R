`RCM.hyperplanes3` <-
function(X,p){
    k <- ncol(X)
    n <- nrow(X)
    N1 <- choose(n,k)
    N2 <- 2^k
    N <- N1*N2
    pN <- floor(p*N)

    samp <- sample(N,pN)
    loc1 <- ceiling(samp/N1) #which block of combinations (1,...,N2)
    loc2 <- samp%%N1+1       #which combination within a block (1,...,N1)
    CombinationsSample <- t(combn(x=n,m=k)[,loc2,drop=FALSE])
   
    binx<-matrix(0,pN,k)
    for(i in 1:k){
      binx[,k-i+1]<-loc1%%2
      loc1<-trunc(loc1/2)
    }
    binx<-2*binx-1
    
    Hyperplanes <- matrix(0, nrow=pN,ncol=k+1)
    for (i in 1:pN){
        Hyperplanes[i,] <- hyperplane(c(binx[i,])*X[CombinationsSample[i,],,drop=FALSE])
    }
    return(Hyperplanes)
}
