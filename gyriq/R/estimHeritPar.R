estimHeritPar <- function(rawResid, Delta, Phi, blkID)
## Finds the maximum likelihood estimator of the polygenic heritability 
## parameter. 'rawResid' is a vector of raw residuals, 'Delta' a censoring 
## indicator, 'Phi' a kinship matrix, and 'blkID' a vector identifying 
## correlated groups of observations.
{
    herit <- optimize(likeHeritPar, interval=c(0, 1), rawResid=rawResid, 
                      Delta=Delta, Phi=Phi, blkID=blkID)$min
    return(herit)
}