calcRegCoeffs <-
function(outList, hBase=1, thin=1, M0=outList$Mcmc$M0, grLabels=paste("Group", 1:outList$Prior$H), printHPD=TRUE, plotPaths=TRUE, plotACFs=TRUE) {

    mnlRes <- matrix(0, dim(outList$Beta.m)[1], (outList$Prior$H-1))

    myc <- 1
    
    myRegCoeffsList <- list()

    for ( h in setdiff(1:outList$Prior$H, hBase) ) {  
  
        tempTable <- rbind( postMean=apply(outList$Beta.m[,h,seq(M0, outList$Mcmc$M, thin)]-outList$Beta.m[,hBase,seq(M0, outList$Mcmc$M, thin)], 1, mean),
                            postSd=apply(outList$Beta.m[,h,seq(M0, outList$Mcmc$M, thin)]-outList$Beta.m[,hBase,seq(M0, outList$Mcmc$M, thin)], 1, sd), 
                            apply(outList$Beta.m[,h,seq(M0, outList$Mcmc$M, thin)]-outList$Beta.m[,hBase,seq(M0, outList$Mcmc$M, thin)], 1, boa.hpd, alpha=0.05) )   
                                
        rownames(tempTable) <- c("Post Exp", "Post Sd", "HPD Lower B", "HPD Upper B") 
        
        myRegCoeffsList[[h]] <- tempTable
    
        if ( printHPD ) {    
            print(xtable(t(tempTable), digits=c(1, rep(5, 4)),
                    caption=paste("Posterior expectation, standard deviation and HPD interval for the MNL regression coefficients: effect of a 
                        certain covariate on the log odds ratio of belonging to group", h, "(", grLabels[h],
                        ") instead of belonging to (baseline) group", hBase, "(", grLabels[hBase], ") (HPD: alpha=0.05)."),
                    label=paste("tab:regCoeffsHPD", h, sep=""))) 
            cat("\n")    
        }
    
        mnlRes[, myc] <- matrix( paste( round((t(tempTable))[,1],5)," (", round((t(tempTable))[,2],5), ")", sep=""), dim(outList$Beta.m)[1], 1) 
     
        myc <- myc + 1
    }

    colnames(mnlRes) <- grLabels[-hBase]
    rownames(mnlRes) <- colnames(tempTable)
    
    myRegCoeffsList$regCoeffsAll <- mnlRes

    print( xtable( mnlRes, digits=c(1, rep(4, outList$Prior$H-1)),
            caption=paste("Posterior expectation and standard deviation of the MNL regression coefficients: effect of a 
                           certain covariate on the log odds ratio of belonging to the indicated group 
                           instead of belonging to (baseline) group", hBase, "(", grLabels[hBase], ")."),
            label="tab:regCoeffs" ) )
            
    if ( plotPaths ) {
    
        dev.new(width=10,height=7) 
        par(mfrow=c(outList$Prior$H - 1, 1), mai=c(0.4, 0.9, 0.1, 0.1), omi=c(0.1,0.0,0.1,0.0)) # c(bottom, left, top, right)  
        for ( h in 2:outList$Prior$H ) {
            matplot(t(outList$Beta.m[,h,]), type="l", lty=1, col=1:8, ylab=grLabels[h], xlab=paste(dimnames(outList$Data$X)[[2]], collapse=", "))   
        }

    }
    
    if ( plotACFs ) {
    
        dev.new(width=7,height=10) 
        par(mfrow=c(min(22,ncol(outList$Data$X)),outList$Prior$H-1), mai=c(0.0, 0.3, 0.1, 0.1), omi=c(0.05, 0.2, 0.25, 0.0))  # c(bottom, left, top, right)
        for ( xv in 1:ncol(outList$Data$X) ) {
            for (h in 2:outList$Prior$H) { 
                acf( outList$Beta.m[xv, h, seq(outList$Mcmc$M0, outList$Mcmc$M, thin)], xaxt="n")
                if ( xv==1 ) mtext(grLabels[h], side = 3, line=1, cex= 1.2 )
                if ( h==2 ) mtext(colnames(outList$Data$X)[xv], side = 2, line=2, cex=if (ncol(outList$Data$X) <= 5 ) 0.8 else 0.5)
            }
        }
 
    }   
    
    return( invisible( myRegCoeffsList ) )

}
