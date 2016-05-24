DHGLMMODELING <-
function(Model="mean",Link=NULL,LinPred="constant",RandDist=NULL,
Offset=NULL,LMatrix=NULL,LinkRandVariance=NULL,LinPredRandVariance=NULL,
RandDistRandVariance="gaussian",LinkRandVariance2=NULL,LinPredRandVariance2=NULL) {
    if (Model=="mean" && is.null(Link)) Link="identity"
    if (Model=="dispersion" && is.null(Link)) Link="log"
    res<-list(Model,Link,LinPred,RandDist,Offset,LMatrix,LinkRandVariance,LinPredRandVariance,
              RandDistRandVariance,LinkRandVariance2,LinPredRandVariance2)
    return(res)
}
