calcTransProbs <-
function(outList, estGroupSize, thin=1, M0=outList$Mcmc$M0, grLabels=paste("Group", 1:outList$Prior$H), printXtable=FALSE, printSd=FALSE, printTogether=TRUE, plotPaths=TRUE, plotPathsForE=TRUE) {

    stopifnot( length(estGroupSize) == outList$Prior$H )
    
    existsXim <- is.element("xi.m", names(outList))

    estTransProb <- 
        if (outList$Prior$H > 1) {        
             if ( existsXim ) { apply( outList$xi.m[seq(M0, outList$Mcmc$M, thin),,,],   c(2, 3, 4), mean) } else { apply( outList$xi_h_m[,,,seq(M0, outList$Mcmc$M, thin)], c(1, 2, 3), mean) }
        } else { 
            if ( existsXim ) { 
                array( apply( outList$xi.m[seq(M0, outList$Mcmc$M, thin),,,],  c(2, 3), mean), c(outList$K+1,outList$K+1,outList$Prior$H)) 
            } else { 
                array( apply(outList$xi_h_m[,,,seq(M0, outList$Mcmc$M, thin)], c(1, 2), mean), c(outList$K+1,outList$K+1,outList$Prior$H))  
            }
        }
        
    dimnames(estTransProb) <- dimnames(outList$Data$dataFile)
    dimnames(estTransProb)[[3]] <- grLabels
    
    estTransProbSd <- 
        if (outList$Prior$H > 1) {
            if ( existsXim ) { apply( outList$xi.m[seq(M0, outList$Mcmc$M, thin),,,], c(2, 3, 4), sd)  } else { apply(outList$xi_h_m[,,,seq(M0, outList$Mcmc$M, thin)], c(1, 2, 3), sd) }
        } else { 
            if ( existsXim ) { 
                array( apply( outList$xi.m[seq(M0, outList$Mcmc$M, thin),,,],  c(2, 3), sd), c(outList$K+1,outList$K+1,outList$Prior$H)) 
            } else {  
                array( apply(outList$xi_h_m[,,,seq(M0, outList$Mcmc$M, thin)], c(1, 2), sd), c(outList$K+1,outList$K+1,outList$Prior$H))
            }
        }        
        
    dimnames(estTransProbSd) <- dimnames(outList$Data$dataFile)
    dimnames(estTransProbSd)[[3]] <- grLabels
     
    if ( printXtable ) {
        for (h in 1:outList$Prior$H) {
            print( xtable( estTransProb[,,h], digits=rep(5, outList$K+2), 
                caption=paste("Posterior expectation of the average transition matrix of group", h, "(", grLabels[h], ") (rel. groupsize: ", round(estGroupSize[h], 4), ")" ), 
                label=paste("tab:trMat", h, sep="")))
            cat("\n")
        }
    }
    
    if ( printSd ) {
        for (h in 1:outList$Prior$H) { 
            print( xtable( 100*estTransProbSd[,,h], digits=rep(5, outList$K+2), 
                caption=paste("Posterior standard deviations (multiplied by 100) of the average transition matrix of group", h, "(", grLabels[h], ")" ), 
                label=paste("tab:trMatSd", h, sep="")))
            cat("\n")
        }
    }    
    
    if ( printTogether ) {              
        for ( h in 1:outList$Prior$H) {         
            matrix1 <- round(estTransProb[,,h], 4)
            matrix2 <- round(100*estTransProbSd[,,h], 3)
            print(xtable( matrix( paste( matrix1," (", matrix2, ")", sep=""), outList$K+1, outList$K+1, dimnames=dimnames(outList$Data$dataFile)[1:2] ), 
                caption=paste("Posterior expectation and posterior standard deviations (multiplied by 100) of the average transition matrix of group", h, "(", grLabels[h], ")" ), 
                label=paste("tab:trMatAndSd", h, sep="")))
            cat("\n")
        }        
    }

    if ( plotPaths ) {    
        upperBoundForXi <- 1.0
        dev.new(width=10,height=7)
        par(mfrow=c(outList$K+1, outList$K+1), mai=c(0.1, 0.1, 0.1, 0.1), omi=c(0.2,0.2,0,0.0)) # c(bottom, left, top, right)
        ylim <- c(0, upperBoundForXi)
        for ( j in 1:(outList$K+1) ) {
            for ( k in 1:(outList$K+1) ) { # draw paths from 1 to M, each group in another color, thinning parameter used        
                plot( if ( existsXim ) { outList$xi.m[seq(1, outList$Mcmc$M, thin), j, k, 1] } else { outList$xi_h_m[j, k, 1, seq(1, outList$Mcmc$M, thin)] },   
                        type="l", col=2, ylim=ylim, xlab="", ylab="", main=NULL, xaxt="n", yaxt=if (k!=1) "n" ) 
                if (j==outList$K+1) axis(1, at=c(0, outList$Mcmc$M0,outList$Mcmc$M/2,outList$Mcmc$M)/thin, labels=paste(c(0, outList$Mcmc$M0, outList$Mcmc$M/2, outList$Mcmc$M)))
                if (outList$Prior$H > 1) for (h in 2:outList$Prior$H) lines( if ( existsXim ) { outList$xi.m[seq(1, outList$Mcmc$M, thin), j, k, h ] } else { outList$xi_h_m[j, k, h, seq(1, outList$Mcmc$M, thin)] }, 
                        type="l", col=h + 1)
            }
        }    
    }    
    
    if ( is.element("e_h_m", names(outList)) & plotPathsForE ) {    
        lowBoundForE <- 45
        dev.new(width=10,height=7) 
        par(mfrow=c(outList$K+1, outList$K+1), mai=c(0.1, 0.1, 0.1, 0.1), omi=c(0.2,0.2,0,0.0)) # c(bottom, left, top, right)
        for ( j in 1:(outList$K+1) ) {
            for ( k in 1:(outList$K+1) ) { # scale of y-axis at least up to lowBoundForE (to be specified)
                ylim <- c(0, max(lowBoundForE, max(outList$e_h_m) + 3))
                plot( outList$e_h_m[j, k, 1, seq(1, outList$Mcmc$M, thin)], type="l", col=2, ylim=ylim, xlab="", ylab="", main=NULL, xaxt="n", yaxt=if (k!=1) "n" )
                if (j==outList$K+1) axis(1, at=c(0, outList$Mcmc$M0,outList$Mcmc$M/2,outList$Mcmc$M)/thin, labels=paste(c(0, outList$Mcmc$M0, outList$Mcmc$M/2, outList$Mcmc$M)))
                if (outList$Prior$H > 1) for (h in 2:outList$Prior$H) lines( outList$e_h_m[j, k, h, seq(1, outList$Mcmc$M, thin)], type="l", col=h + 1)
            }
        }
    }
        
    return( invisible( list( estTransProb=estTransProb, estTransProbSd=estTransProbSd ) ) )
}
