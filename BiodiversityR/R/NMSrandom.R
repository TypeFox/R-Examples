`NMSrandom` <-
function(x,perm=100,k=2,stressresult=F,method="isoMDS"){
#    if (!require(MASS)) {stop("Requires package MASS")}
    minstress <- 100
    stress <- array(dim=perm)
    for (j in 1:perm) {
        if (method=="isoMDS") {result <- MASS::isoMDS(x,initMDS(x,k=k),k=k,maxit=1000,trace=F)}
        if (method=="sammon") {result <- MASS::sammon(x,initMDS(x,k=k),k=k,niter=1000,trace=F)}
        stress[j] <- result$stress
        if (result$stress < minstress) {
            minstress <- result$stress
            minresult <- result
        }
    }
    rownames(minresult$points) <- rownames(as.matrix(x))
    if (stressresult==F) {return(minresult)}
    if (stressresult==T) {return(stress)}
}

