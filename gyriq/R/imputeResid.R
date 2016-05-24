imputeResid <- function(rawResid, Delta, Gam, m, blkID)
## Imputes 'm' vectors of residuals for censored failure times and returns the
## resulting completed vector of residuals 'compResid' once standardized.
##
## 'rawResid': vector of raw residuals
## 'Delta': vector containing the censoring indicator
## 'Gam': estimated covariance matrix of the raw residuals
## 'blkID': vector with entries identifying correlated groups of observations
{
    cn <- estimCondNorm(as.logical(Delta), rawResid, Gam)
        # 'cn' is a list from which the mean and variance of the conditional
        # normal distribution of the raw residuals of the censored failure times
        # are retrieved.
    censC <- Delta == 0
    censTot <- sum(censC)
    censInd <- 1:censTot
    
    blkIDcens <- blkID[censC]
    blkUniqCens <- unique(blkIDcens)
    blkTot <- length(blkUniqCens)
    
    impList <- apply(as.matrix(1:blkTot), MARGIN=1, FUN=imputeBlkResid, 
                     blkUniqCens=blkUniqCens, blkIDcens=blkIDcens, cn=cn, m=m)
    impVec <- unlist(impList)
    impResid <- matrix(impVec, nrow=censTot, ncol=m, byrow=TRUE)
    
    compResid <- numeric(length(rawResid))
    compResid[cn$obsInd] <- cn$obsResid
    compResid[cn$censInd] <- apply(impResid, 1, mean)
    
    GamEig <- eigen(Gam)
    GamEigVal <- GamEig$values
    GamSqrt <- GamEig$vectors %*% diag(sqrt(GamEigVal)) %*% 
        solve(GamEig$vectors)
    rho <- mean(compResid^2)
        ## Scale parameter which reflects the fact that we are using multiple
        ## imputed values rather than real observations.
    return(as.vector(solve(GamSqrt) %*% compResid / sqrt(rho)))
}