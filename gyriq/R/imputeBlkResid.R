imputeBlkResid <- function(blkIndex, blkUniqCens, blkIDcens, cn, m)
## Imputes 'm' vectors of residuals for the censored failure times of a 
## correlated group of observations (identified by 'blkIndex', 'blkUniqCens' and
## 'blkIDcens').
## 
## 'cn' is a list from which the raw residuals of the censored failure times are
## retrieved along with the mean and variance of their conditional normal 
## distribution.
{
    sbjInd <- blkIDcens == blkUniqCens[blkIndex]
    sbjMean <- cn$meanCond[sbjInd]
    sbjVar <- cn$varCond[sbjInd, sbjInd]
    sbjLower <- cn$censResid[sbjInd]
    
    if (length(sbjMean) == 1) {
        sbjUnif <- 1 - pnorm(rnorm(m, mean=sbjMean, sd=sqrt(sbjVar)))
    } else {
        sbjUnif <- 1 - pnorm(mvtnorm::rmvnorm(m, mean=sbjMean, sigma=sbjVar))
      }
    sbjLowerP <- matrix(rep(1 - pnorm(sbjLower), m), nrow=m,
                        ncol=length(sbjLower), byrow=TRUE)
    sbjResid <- qnorm(1 - sbjLowerP * sbjUnif)
        # Generates the vectors of imputed values with the restriction to be
        # larger than the original censored values, componentwise.
    return(as.matrix(sbjResid))
}