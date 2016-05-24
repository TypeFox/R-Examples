estimCondNorm <- function(blkDeltaC, blkResid, Gam)
## Computes the mean and the variance of a conditional normal distribution used
## when the log-likelihood function of the raw residuals is expressed as a 
## density function times a conditional survival function.
##
## 'blkResid' and 'blkDeltaC' are respectively the raw residuals and the 
## censoring indicator of the group of correlated observations for which the 
## mean and variance are computed. 'Gam' is the variance of 'blkResid'.
{
    obsInd <- seq(along.with=blkDeltaC)[blkDeltaC]
    censInd <- seq(along.with=blkDeltaC)[!blkDeltaC]
    obsResid <- blkResid[obsInd]
    censResid <- blkResid[censInd]
 
    Gam11 <- Gam[obsInd, obsInd]
    Gam10 <- Gam[obsInd, censInd]
    Gam01 <- Gam[censInd, obsInd]
    Gam00 <- Gam[censInd, censInd]
    
    invm <- Gam01 %*% solve(Gam11)
    meanCond <- as.vector(invm %*% obsResid)
    varCond <- as.matrix(Gam00 - invm %*% Gam10)
    
    return(list(meanCond=meanCond, varCond=varCond, obsInd=obsInd, 
                censInd=censInd, obsResid=obsResid, censResid=censResid, 
                Gam11=Gam11))
}