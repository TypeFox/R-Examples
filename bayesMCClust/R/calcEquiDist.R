calcEquiDist <-
function(outList, thin=1, maxi=50, M0=outList$Mcmc$M0, grLabels=paste("Group", 1:outList$Prior$H), printEquiDist=TRUE, plotEquiDist=TRUE) {

    m_mi <- function(mi, maxim=maxi, M=outList$Mcmc$M, M0=outList$Mcmc$M0) {
        stopifnot( M - (maxim - 1)*thin >= M0 )
        return( M - (maxim - mi)*thin)     
    }

    equiDist <- function(K, trMat) {
        A <- cbind(diag(1,K+1,K+1) - trMat, 1)
        (t(A) %*% solve( A %*% t(A)))[K+2,]
    }

    edArray <- array(0, c( outList$K+1, outList$Prior$H, maxi ))

    for (mi in 1:maxi) edArray[ , , mi ] <- sapply(1:outList$Prior$H, function(h) equiDist(K=outList$K, if ( is.element("xi.m", names(outList)) ) { outList$xi.m[m_mi(mi),,,h] } else { outList$xi_h_m[,,h,m_mi(mi)] } ))   

    ed <- apply(edArray, c(1,2), mean)
    rownames(ed) <- paste(as.numeric(dimnames(outList$Data$dataFile)[[1]]), ": ",sep="")
    colnames(ed) <- grLabels     
    
    if ( printEquiDist ) {
        print(xtable(ed, digits=rep(5,outList$Prior$H+1), 
            caption=paste("Steady state (stationary distribution) of the transition matrix in the various clusters (using only last", maxi, "of each", thin,"-th draw)."), 
            label="tab:ed"))
        cat("\n")
    } 
    
    if ( plotEquiDist ) {
    
        dev.new()
        par(mai=c(0.5, 0.6, 0.5, 0.1), omi=c(0.2, 0.2, 0.1, 0.1))
        barplot2(ed, beside=FALSE, col = rainbow(outList$K+1), cex.names=1.3, las=1, cex.axis=1.0, main="Stationary distribution") 

    }   
    
    return( invisible( ed ) )

}
