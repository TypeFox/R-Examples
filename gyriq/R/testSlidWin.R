testSlidWin <- function(indexSW, compResid, G, w, ker, asv, method, starResid, 
                        bsw, tsw, pos, fileOut)
## Calculates the score statistic and p-value of the kinship-adjusted survival
## SNP-set association test for one sliding window identified by the index
## 'indexSW' among the 'length(bsw)' possible sliding windows.
##
## 'compResid': completed vector of residuals
## 'G': matrix containing a set of SNPs
## 'w': vector of weights for the SNPs
## 'asv': (default=NULL) number of approximate eigenvalues to be estimated
##      for the kernel matrix if the implicitly-restarted Lanczos 
##      bidiagonalization is used
## 'method': (default="davies") procedure used to obtain the p-value
## 'starResid': (default=NULL) matrix of permuted residuals used to obtain the
##      p-value if a permutation procedure is employed
## 'bsw': lower bound of the sliding window
## 'tsw': upper bound of the sliding window
## 'pos': vector of SNP positions (used for the output only)
## 'fileOut': a string containing the name and path of the output file where the
##      results are to be printed
{
    indSNP <- bsw[indexSW]:tsw[indexSW]
    nsw <- length(indSNP)
    gyr <- try(testGyriq(compResid=compResid, G=G[, indSNP], w=w[indSNP], 
                         ker=ker, asv=asv, method=method, starResid=starResid), 
               silent=TRUE)
    if (inherits(gyr, "try-error")) {
        gyrMes <- attr(x=gyr, which="condition")$message
        gyr <- data.frame(score=-99, pVal=-99)
    } else {
        gyrMes <- "OK"
    }
    resSW <- c(indSNP[1], indSNP[nsw], nsw, pos[indSNP[1]], pos[indSNP[nsw]], 
        signif(gyr$score), signif(gyr$pVal), gyrMes)
    write.table(as.data.frame(t(resSW)), file=fileOut, append=TRUE, quote=FALSE,
        sep="\t", dec=".", row.names=FALSE, col.names=FALSE)
    print(paste("Processing sliding window ", indexSW, sep=""))
    return(resSW)
}