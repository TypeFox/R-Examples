`RCM.hyperplanes2` <-
function(X,p){
    X2 <- rbind(X,-X)
    k <- ncol(X2)
    n <- nrow(X2)
    N <- choose(n,k)
    pN <- floor(p*N)
    CombinationsSample <- createHyperplaneSample(N,pN,n,k)
    Hyperplanes <- matrix(0, nrow=pN,ncol=k+1)
    for (i in 1:pN){
        Hyperplanes[i,] <- hyperplane(X2[CombinationsSample[i,],,drop=FALSE])
    }
    return(Hyperplanes)
}
