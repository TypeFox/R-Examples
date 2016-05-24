calcAllocationsMCCExt <-
function(outList, thin=1, maxi=50, M0=outList$Mcmc$M0) {

    m_mi <- function(mi, maxim=maxi, M=outList$Mcmc$M, M0=outList$Mcmc$M0) {
        stopifnot( M - (maxim - 1)*thin >= M0 )
        return( M - (maxim - mi)*thin)     
    }

    sProbsNormSum <- matrix(0, outList$N, outList$Prior$H)
    sProbsGrSize <- matrix(0, maxi, outList$Prior$H) 
    
    for (mi in 1:maxi) { 
        XBetak <- crossprod(t(outList$Data$X), outList$Beta.m[,,m_mi(mi)])
        logit.temp <- exp( XBetak ) / rowSums( exp( XBetak ) ) 
        sProbs <- 
        matrix(mapply(function(i,h) prod(outList$xi.m[m_mi(mi),,,h]^outList$Njk.i[,,i]),rep(1:outList$N,each=outList$Prior$H),1:outList$Prior$H ),outList$N,outList$Prior$H,byrow=TRUE)*logit.temp
        sProbsNorm  <-  sProbs/rowSums(sProbs)    
        sProbsNormSum <- sProbsNormSum + sProbsNorm    
        sProbsGrSize[mi,] <- colMeans(sProbsNorm)    
    
        if ( identical(all.equal(mi %% 10, 0), TRUE) ) {
            cat("", mi, "of", maxi, "iterations done...\n")
            flush.console()
        }
    
    }
    
    estGrSizes <- colMeans(sProbsGrSize)
    estGroupSize <- estGrSizes  
    class <- max.col(sProbsNormSum) 
    cpTemp <- sProbsNormSum 
    cp2 <- sProbsNormSum/maxi

    allocList <- list(estGroupSize=estGroupSize, class=class, classProbs=cp2)
    
    return( invisible( allocList ) )

}
