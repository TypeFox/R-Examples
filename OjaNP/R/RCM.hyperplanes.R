`RCM.hyperplanes` <-
function(X,p){
    k <- ncol(X)
    n <- nrow(X)
    N <- choose(n,k)
    pN <- floor(p*N)
    CombinationsSample <- createHyperplaneSample(N,pN,n,k)
    Hyperplanes <- matrix(0, nrow=pN,ncol=k+1)
    for (i in 1:pN){
        Hyperplanes[i,] <- hyperplane(X[CombinationsSample[i,],,drop=FALSE])
    }
    return(Hyperplanes)
}
