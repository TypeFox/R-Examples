likeHeritPar <- function(herit, rawResid, Delta, Phi, blkID)
## Computes the value of the log-likelihood function of the raw residuals 
## 'rawResid' when the polygenic heritability parameter is equal to 'herit'.
##
## The censoring indicator 'Delta' and the kinship matrix 'Phi' enters in the
## definition of the log-likelihood function. 'blkID' is a vector with entries
## identifying correlated groups of observations.
{
    blkUniqVec <- unique(blkID)
    blkTot <- length(blkUniqVec)
    likeList <- apply(as.matrix(1:blkTot), MARGIN=1, FUN=calcLikeBlk, 
                      blkUniqVec=blkUniqVec, herit=herit, rawResid=rawResid,
                      Delta=Delta, Phi=Phi, blkID=blkID)
    likeVec <- unlist(likeList)
    likeVal <- - sum(likeVec)
    return(likeVal)
}