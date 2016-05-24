calcAllocationsMCC <-
function(outList, thin=1, maxi=50, M0=outList$Mcmc$M0, plotPathsForEta=TRUE) {

    m_mi <- function(mi, maxim=maxi, M=outList$Mcmc$M, M0=outList$Mcmc$M0) {
        stopifnot( M - (maxim - 1)*thin >= M0 )
        return( M - (maxim - mi)*thin)     
    }

    # calculation of classification probabilities
    sProbsNormSum <- matrix(0, outList$N, outList$Prior$H)
    for (mi in 1:maxi) { 
        sProbs <- 
            matrix(mapply(function(i,h) prod(outList$xi.m[m_mi(mi),,,h]^outList$Njk.i[,,i]), rep(1:outList$N, each=outList$Prior$H), 1:outList$Prior$H ), outList$N, outList$Prior$H,
                byrow=TRUE)*matrix(outList$eta.m[,m_mi(mi)], outList$N, outList$Prior$H, byrow=TRUE)
        sProbsNorm  <-  sProbs/rowSums(sProbs)    
        sProbsNormSum <- sProbsNormSum + sProbsNorm
    
        if ( identical(all.equal(mi %% 10, 0), TRUE) ) {
            cat("", mi, "of", maxi, "iterations done...\n")
            flush.console()
        }

    
    }
    class <- max.col(sProbsNormSum)  
    cpTemp <- sProbsNormSum 
    cp2 <- sProbsNormSum/maxi 

    estGroupSize <- 
        if (outList$Prior$H > 1) { 
            apply( outList$eta.m[ , seq(M0, outList$Mcmc$M, thin)], 1, mean ) 
        } else { 
            mean( outList$eta.m[ , seq(M0, outList$Mcmc$M, thin)] )
        }

    if ( plotPathsForEta ) {    
        dev.new(width=10,height=7) 
        par(mfrow=c(1, 1))
        if (outList$Prior$H > 1) matplot( t(outList$eta.m[, seq(1, outList$Mcmc$M, thin)]), type="l", col=2:(outList$Prior$H+1), xlab="", ylab="", lty=1, lwd=0.5, 
            main="MCMC draws for group sizes", xaxt="n", ylim=c(0, min(1, max(outList$eta.m)+0.05))) else matplot( outList$eta.m[, seq(1, outList$Mcmc$M, thin)], 
            type="l", col=2:(outList$Prior$H+1), xlab="", ylab="", lty=1, lwd=0.5, main="MCMC draws for group sizes", xaxt="n", ylim=c(0, min(1, max(outList$eta.m)+0.05)))
        axis(1, at=c(0, M0,outList$Mcmc$M/2,outList$Mcmc$M)/thin, labels=paste(c(0, M0, outList$Mcmc$M/2, outList$Mcmc$M)))        
    
    }

    allocList <- list(estGroupSize=estGroupSize, class=class, classProbs=cp2)
    
    return( invisible( allocList ) )

}
