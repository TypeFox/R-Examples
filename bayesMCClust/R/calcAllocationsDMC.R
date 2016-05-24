calcAllocationsDMC <-
function(outList, thin=1, maxi=50, M0=outList$Mcmc$M0, plotPathsForEta=TRUE) {

    m_mi <- function(mi, maxim=maxi, M=outList$Mcmc$M, M0=outList$Mcmc$M0) {
        stopifnot( M - (maxim - 1)*thin >= M0 )
        return( M - (maxim - mi)*thin)     
    }

    margDataLL_m <- function(i, h, m, Njk.i=outList$Njk.i, e_h_m=outList$e_h_m, AA=AA, BB=BB, plotPathsForEta=TRUE) {
        CC <- sum( lgamma( Ne <- (Njk.i[,,i] + e_h_m[,,h,m])   ) )
        DD <- sum( lgamma( rowSums( Ne ) ) )
        return( exp(AA[h] - BB[h] + CC - DD) )
    }

    # calculation of classification probabilities
    margDataLikeli <- matrix(0, outList$N, outList$Prior$H)
    sProbsNormSum <- matrix(0, outList$N, outList$Prior$H)
    if (outList$Prior$H > 1) {
        for (mi in 1:maxi) {  
            AA <- apply( lgamma(apply(outList$e_h_m[,,,m_mi(mi)], 3, rowSums)), 2, sum )  
            BB <- apply( lgamma(outList$e_h_m[,,,m_mi(mi)]), 3, sum )
            for ( i in 1:outList$N) {
                for ( h in 1:outList$Prior$H ) {
                    margDataLikeli[i, h] <- margDataLL_m(i=i, m=m_mi(mi), h=h, Njk.i=outList$Njk.i, e_h_m=outList$e_h_m, AA=AA, BB=BB)
                }
            }  
            sProbs  <- margDataLikeli*matrix(outList$eta_m[m_mi(mi),], outList$N, outList$Prior$H, byrow=TRUE) 

            sProbsNorm  <-  sProbs/rowSums(sProbs)    
            sProbsNormSum <- sProbsNormSum + sProbsNorm    
    
            if ( identical(all.equal(mi %% 10, 0), TRUE) ) {
                cat("", mi, "of", maxi, "iterations done...\n")
                flush.console()
            }

        }
    }
    class <- max.col(sProbsNormSum)
    cpTemp <- sProbsNormSum 
    cp2 <- if (outList$Prior$H > 1) { sProbsNormSum/maxi } else { matrix(1, outList$N, outList$Prior$H) }

    estGroupSize <- 
        if (outList$Prior$H > 1) { 
            apply( outList$eta_m[seq(M0, outList$Mcmc$M, thin),], 2, mean ) 
        } else { 
            mean( outList$eta_m[seq(M0, outList$Mcmc$M, thin),] ) 
        }
        
    if ( plotPathsForEta ) {    
        dev.new(width=10,height=7) 
        par(mfrow=c(1, 1))
        if (outList$Prior$H > 1) matplot( (outList$eta_m[seq(1, outList$Mcmc$M, thin), ]), type="l", col=2:(outList$Prior$H+1), xlab="", ylab="", lty=1, lwd=0.5, 
            main="MCMC draws for group sizes", xaxt="n", ylim=c(0, min(1, max(outList$eta_m)+0.05))) else matplot( outList$eta_m[ seq(1, outList$Mcmc$M, thin),], 
            type="l", col=2:(outList$Prior$H+1), xlab="", ylab="", lty=1, lwd=0.5, main="MCMC draws for group sizes", xaxt="n", ylim=c(0, min(1, max(outList$eta_m)+0.05)))
        axis(1, at=c(0, M0,outList$Mcmc$M/2,outList$Mcmc$M)/thin, labels=paste(c(0, M0, outList$Mcmc$M/2, outList$Mcmc$M)))    
    }

    allocList <- list(estGroupSize=estGroupSize, class=class, classProbs=cp2)
    
    return( invisible( allocList ) )

}
