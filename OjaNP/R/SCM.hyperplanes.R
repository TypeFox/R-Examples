`SCM.hyperplanes` <-
function(X,center,p,...){
    # center has to be 'numeric'
    k <- ncol(X)
    n <- nrow(X)    
    N <- choose(n,k-1)
    pN <- floor(p*N)
    CombinationsSample <- createHyperplaneSample(N,pN,n,k-1)
    Hyperplanes <- matrix(0, nrow=pN, ncol=k+1)
    for (i in 1:pN){
          Hyperplanes[i,] <- hyperplane(rbind(center,X[CombinationsSample[i,],,drop=FALSE])) 
    }
    return(Hyperplanes)
}
