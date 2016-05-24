calcVariationDMC <-
function(outList, thin=1, maxi=50, M0=outList$Mcmc$M0, grLabels=paste("Group", 1:outList$Prior$H), printVarE=FALSE, printUnobsHet=FALSE, printUnobsHetSd=FALSE, 
                             printUnobsHetAll=FALSE, printAllTogether=TRUE) {

    m_mi <- function(mi, maxim=maxi, M=outList$Mcmc$M, M0=outList$Mcmc$M0) {
        stopifnot( M - (maxim - 1)*thin >= M0 )
        return( M - (maxim - mi)*thin)     
    }

    # elementwise variances
    var_e_Array <- array( 0, c(outList$K+1, outList$K+1, outList$Prior$H, maxi) )
    for (mi in 1:maxi) {
        for (h in 1:outList$Prior$H) {
            var_e_Array[,,h,mi] <-  outList$e_h_m[,,h,m_mi(mi)]*( rowSums(outList$e_h_m[,,h,m_mi(mi)])-outList$e_h_m[,,h,m_mi(mi)]) /
                                                        ((rowSums(outList$e_h_m[,,h,m_mi(mi)])^2)*(1+rowSums(outList$e_h_m[,,h,m_mi(mi)])))
        }
    }
    var_e <- apply(var_e_Array, c(1,2,3), mean)
    
    dimnames(var_e) <- dimnames(outList$Data$dataFile)
    dimnames(var_e)[[3]] <- grLabels

    if ( printVarE ) {
        for (h in 1:outList$Prior$H) { 
            print( xtable(var_e[,,h]*10^(4), digits=rep(3, outList$K+2), 
                caption=paste("Posterior expectation of the variance of the individual transition probabilities (in percent) in group", h, "(", grLabels[h], 
                              ") (using only last", maxi, "of each", thin,"-th draw)."),
                label=paste("tab:varIndTransProb",h,sep=""))) 
            cat("\n")
        }
    }

    # amount of heterogeneity
    het_Array <-  array(0, c( outList$K+1, outList$Prior$H, maxi ))
    for (mi in 1:maxi) het_Array[ , , mi ] <- (1+apply(array(outList$e_h_m[,,,m_mi(mi)], c(outList$K+1, outList$K+1, outList$Prior$H)), c( 1, 3), sum))^(-1)
    het <- apply(het_Array, c(1,2), mean)
    dimnames(het) <- list(dimnames(outList$Data$dataFile)[[1]], grLabels) 
    
    if ( printUnobsHet ) {
        print(xtable(cbind(100*het), digits=rep(5,outList$Prior$H+1), 
            caption=paste("Posterior expectation of the row-specific unobserved heterogeneity measure in each group and multiplied by 100."), 
            label="tab:unobsHet"))
        cat("\n")
    }

    #==============================================================
    hetsd <- apply(het_Array, c(1,2), sd)
    dimnames(hetsd) <- list(dimnames(outList$Data$dataFile)[[1]], grLabels) 
    
    if ( printUnobsHetSd ) {
        print(xtable(cbind(100*hetsd), digits=rep(5,outList$Prior$H+1), 
            caption=paste("Posterior standard deviation of the row-specific unobserved heterogeneity measure in each group multiplied by 100 (using only last", maxi, "of each", thin,"-th draw)."), 
            label="tab:unobsHetSd"))
        cat("\n")
    }

    # =============================================================
    matrix1 <- round(cbind(100*het), 3) 
    matrix2 <- round(cbind(100*hetsd), 3)
    myTempMat <- matrix( paste( matrix1," (", matrix2, ")", sep=""), outList$K+1, outList$Prior$H, dimnames=list( dimnames(outList$Data$dataFile)[[1]], grLabels))
    
    if ( printUnobsHetAll ) {    
        print( xtable( myTempMat,
            caption=paste("Posterior expectation and, in parenthesis, posterior standard deviation of the row-specific unobserved heterogeneity measure in each group, each multiplied by 100 
                           (using only last", maxi, "of each", thin,"-th draw)."),
            label="tab:unobsHetAll"))
        cat("\n")
    }

    # =============================================================
    if ( printAllTogether ) {    
        for (h in 1:outList$Prior$H) {
            print( xtable( cbind( round( var_e[,,h]*10^(4), 3), unobsHet=myTempMat[,h]) , digits=rep(4, outList$K+2+1), 
                caption=paste("Posterior expectation of the variance of the individual transition probabilities (in percent) in group", h, "(", grLabels[h], ") and in the last column: 
                    posterior expectation and, in parenthesis, posterior standard deviation of the row-specific unobserved heterogeneity measure in each group, each multiplied by 100 
                    (using only last", maxi, "of each", thin,"-th draw)."),
                label=paste("tab:varIndTransProbANDunobsHet", h, sep="")))
            cat("\n")
        }
    }
    
    return( invisible( list( var_e=var_e, het=het, hetsd=hetsd ) ) )

}
