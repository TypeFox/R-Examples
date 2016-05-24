calcLongRunDist <-
function(outList, initialStateData, class, equiDist, thin=1, maxi=50, M0=outList$Mcmc$M0, printLongRunDist=TRUE, grLabels=paste("Group", 1:outList$Prior$H)) {

    "%^%" <- function(mat, pow) { 
        stopifnot(length(pow) == 1, all.equal(pow, round(pow)), nrow(mat) == ncol(mat)) 
        pow <- round(pow) 
        if (pow < 0) { 
            mat <- solve(mat) 
            pow <- abs(pow) 
        } 
        result <- diag(nrow(mat)) 
        while (pow > 0) { 
            result <- result %*% mat 
            pow <- pow - 1 
        } 
        result 
    } 

    m_mi <- function(mi, maxim=maxi, M=outList$Mcmc$M, M0=outList$Mcmc$M0) {
        stopifnot( M - (maxim - 1)*thin >= M0 )
        return( M - (maxim - mi)*thin)     
    }
    
    katLabels <- unique( as.numeric( unlist( dimnames( outList$Data$dataFile)[1:2]), names( table( initialStateData ))))
    
    if ( is.element("xi.m", names(outList)) ) { # mcc
    
        xiArray <- array( 0 , c( outList$K+1, outList$K+1, outList$Prior$H, maxi) )
        for (mi in 1:maxi ) {  # for each thin-th draw from the last ones...
            for ( h in 1:outList$Prior$H ) { # for each group...
                xiArray[,,h,mi] <- outList$xi.m[m_mi(mi),,,h] # mi
            }
        }

        times <- c( 1,2,3,4,5, 10, 15, 20, 25, 30, 50, 100 )  
        xiht <- vector( mode="list", length=max(times) )
        for (t in times) {
            xiTtemp <- apply( xiArray, c(3,4), "%^%", t )   #  to exponentiate 
            xiTmean <- apply( xiTtemp, c(1,2), mean )   # average over draws m    
            xihTemp <- array(0, c(outList$K+1,outList$K+1,outList$Prior$H))    
            for ( h in 1:outList$Prior$H ) {      
                xihTemp[,,h] <- matrix( xiTmean[ , h ], outList$K+1, outList$K+1 ) #  within group    
            }    
            xiht[[t]] <- xihTemp
        }

    } else { # dmc
    

        xiArray <- array( 0 , c( outList$K+1, outList$K+1, outList$N, maxi) )
        for (mi in 1:maxi) { # for each thin-th draw from the last ones...
            for ( i in 1:outList$N ) { # for each individual...
                xi_i <- matrix(0, outList$K+1, outList$K+1)
                for ( jR in 1:(outList$K+1) ) {  #   for each row...
                    xi_i[jR,] <- gtools::rdirichlet(n=1, alpha=outList$e_h_m[jR, , class[i], m_mi(mi) ] + outList$Njk.i[jR,,i] ) # posterior distribution 
                }    
                xiArray[,,i,mi] <- xi_i
            }
        }
    
        times <- c( 1,2,3,4,5, 10, 15, 20, 25, 30, 50, 100 ) 
        xiht <- vector( mode="list", length=max(times) )
        for (t in times) {

            partSize <- min(outList$N, 2000)
      
            loopMin <- 1
            loopMax <- partSize    
            xiTmean <- matrix(0, (outList$K+1)*(outList$K+1), outList$N )
    
            while ( loopMin < outList$N ) {    
        
                xiTtemp <- apply( xiArray[ , , loopMin:loopMax , ], c(3,4), "%^%", t )          
                xiTmean[ , loopMin:loopMax ] <- apply( xiTtemp, c(1,2), mean )        
                rm(xiTtemp)    
                loopMin <- loopMin + partSize
                loopMax <- min( loopMax + partSize, outList$N)    
            }
    
            xihTemp <- array(0, c(outList$K+1,outList$K+1,outList$Prior$H))    
            for ( h in 1:outList$Prior$H ) {      
                xihTemp[,,h] <- matrix( apply( xiTmean[ , class==h ], 1, mean), outList$K+1, outList$K+1) # average within group 
            }    
            xiht[[t]] <- xihTemp
        }
    
    }

    myLongRunDistList <- list()
    
    # long run dist incl Inf
    dev.new(width=10,height=7) 
    par(mfrow=c(outList$Prior$H, 1), mai=c(0.3, 0.6, 0.1, 0.0), omi=c(0.0, 0.0, 0.0, 0.0))  # c(bottom, left, top, right)
    for ( h in 1:outList$Prior$H) {
        timep <- c(0, times) 
        nt <- length(timep)
        tabProgr <- matrix(0, outList$K+1, nt, dimnames = list( paste( katLabels ), c(paste("t=",timep[1:6],sep=""),timep[-(1:6)]) ))
        tabProgr[,1] <- (table(c(katLabels, initialStateData[class==h]))-1)/sum(class==h)
        tabProgr[,2:nt] <- 
            sapply(timep[2:nt], function(t) rbind( (table(c(katLabels, initialStateData[class==h]))-1)/sum(class==h)) %*% ( xiht[[t]][,,h] ))
        tabProgr[,nt] <- equiDist[,h]
        dimnames(tabProgr)[[2]][nt] <- Inf
        barplot2(tabProgr, beside=FALSE, col = rainbow(outList$K+1), ylab=grLabels[h], las=1, cex.main = 1.7, cex.names=1.8, cex.axis=1.4, cex.lab=1.9 ) 
        
        if ( printLongRunDist ) {
        
            print( xtable( tabProgr , digits=rep(4, nt+1), 
                caption=paste("Long-run distribution of group", h, "(", grLabels[h], ")"), 
                label=paste("tab:longRunDist", h, sep="") ) )  
                
            cat("\n")     
        }       
        
        myLongRunDistList[[h]] <- tabProgr
    }
    
    return( invisible( myLongRunDistList ) )   

}
