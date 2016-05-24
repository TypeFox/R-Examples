calcParMatDMC <-
function(outList, thin=1, M0=outList$Mcmc$M0, grLabels=paste("Group", 1:outList$Prior$H), printPar=TRUE) {

    e_h_cP <- 
    if (outList$Prior$H > 1) {
               apply(outList$e_h_m[,,,seq(M0,outList$Mcmc$M,thin)], c(1, 2, 3), mean) 
    } else { 
        array( apply(outList$e_h_m[,,,seq(M0,outList$Mcmc$M,thin)], c(1, 2), mean) , c(outList$K+1,outList$K+1,outList$Prior$H))  
    }
    
    dimnames(e_h_cP) <- dimnames(outList$Data$dataFile)
    dimnames(e_h_cP)[[3]] <- grLabels
    
    if ( printPar ) {    
        for (h in 1:outList$Prior$H) { 
            print( xtable( e_h_cP[,,h], digits=rep(4, outList$K+2), 
                caption=paste("Parameter matrix of group", h, "(", grLabels[h], ")"), 
                label=paste("tab:parMat", h, sep=""))) 
            cat("\n")
        }    
    }
    
    return( invisible( e_h_cP ) )
    
}
