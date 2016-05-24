calcLikeBlk <- function(blkIndex, blkUniqVec, herit, rawResid, Delta, Phi, 
                        blkID)
## Computes the contribution of a correlated group of observations (identified 
## by 'blkIndex', 'blkUniqVec' and 'blkID') to the value of the log-likelihood 
## function of the raw residuals 'rawResid' when the polygenic heritability 
## parameter is equal to 'herit'.
##
## The censoring indicator 'Delta' and the kinship matrix 'Phi' enters in the
## definition of the log-likelihood function.
{
    sbjInd <- blkID == blkUniqVec[blkIndex]
    sbjTot <- sum(sbjInd)
    Imat <- diag(rep(1, sbjTot))
    Gam <- herit * Phi[sbjInd, sbjInd] + (1 - herit) * Imat
    
    blkResid <- rawResid[sbjInd]
    blkDeltaC <- Delta[sbjInd] == 1
    obsTot <- sum(blkDeltaC)
    censTot <- sum(sbjInd) - obsTot
    if (obsTot == 0) {
        contribBlk <- log(mvtnorm::pmvnorm(lower=blkResid, sigma=Gam)[1])
    } else {
        if (censTot == 0) {
            contribBlk <- mvtnorm::dmvnorm(blkResid, sigma=Gam, log=TRUE)
        } else {
            ## Log-likelihood function expressed as a density function times a
            ## survival function of a conditional normal distribution
            cn <- estimCondNorm(blkDeltaC=blkDeltaC, blkResid=blkResid, Gam=Gam)
            contribObs <- mvtnorm::dmvnorm(cn$obsResid, 
                                           sigma=as.matrix(cn$Gam11), log=TRUE)
            contribCens <- log(mvtnorm::pmvnorm(lower=cn$censResid, 
                                                mean=cn$meanCond, 
                                                sigma=cn$varCond)[1])
            contribBlk <- c(contribObs, contribCens)
        }
    }
    return(contribBlk)
}