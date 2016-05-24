calcAllocationsDMCExt <-
function(outList, thin=1, maxi=50, M0=outList$Mcmc$M0) {

    m_mi <- function(mi, maxim=maxi, M=outList$Mcmc$M, M0=outList$Mcmc$M0) {
        stopifnot( M - (maxim - 1)*thin >= M0 )
        return( M - (maxim - mi)*thin)     
    }

    margDataLL_m <- function(i, h, m, Njk.i=outList$Njk.i, e_h_m=outList$e_h_m, AA=AA, BB=BB) {
        CC <- sum( lgamma( Ne <- (Njk.i[,,i] + e_h_m[,,h,m])   ) )
        DD <- sum( lgamma( rowSums( Ne ) ) )
        return( exp(AA[h] - BB[h] + CC - DD) )
    }

    # calculation of classification probabilities
    margDataLikeli <- matrix(0, outList$N, outList$Prior$H)
    sProbsNormSum <- matrix(0, outList$N, outList$Prior$H)
    sProbsGrSize <- matrix(0, maxi, outList$Prior$H)
    for (mi in 1:maxi) {  
        XBetak <- crossprod(t(outList$Data$X), outList$Beta.m[,,m_mi(mi)])   
        logit.temp <- exp( XBetak ) / rowSums( exp( XBetak ) ) 

        AA <- apply( lgamma(apply(outList$e_h_m[,,,m_mi(mi)], 3, rowSums)), 2, sum )  
        BB <- apply( lgamma(outList$e_h_m[,,,m_mi(mi)]), 3, sum )
        for ( i in 1:outList$N) {
            for ( h in 1:outList$Prior$H ) {
                margDataLikeli[i, h] <- margDataLL_m(i=i, m=m_mi(mi), h=h, Njk.i=outList$Njk.i, e_h_m=outList$e_h_m, AA=AA, BB=BB)
            }
        }  
        sProbs  <- margDataLikeli*logit.temp 

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
