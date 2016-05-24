plotTransProbs <-
function(outList, estTransProb, estGroupSize, class, grLabels=paste("Group", 1:outList$Prior$H), plotPooled=TRUE, plotContTable=TRUE, printContTable=TRUE, plotContPooled=TRUE) {
 
    gewTransProb <- estTransProb
    Njk <- apply(outList$Njk.i, c(1, 2), sum)
    relNjk <- Njk/apply(Njk, 1, sum)    
    scal <- 1.0/max(gewTransProb, relNjk) 
    k   <- rep(as.numeric(dimnames( outList$Data$dataFile)[[1]]), outList$K+1)
    j   <- rep(as.numeric(dimnames( outList$Data$dataFile)[[1]]), each = outList$K+1)
    
    dev.new()
    par(mfrow=c(ceiling(outList$Prior$H/2),2), mai=c(0.2, 0.2, 0.2, 0.25))
    sapply(1:outList$Prior$H, function(h) balloonplot( k, j, c(t(gewTransProb[,,h])), 
        scale.method = c("volume"), main=paste(grLabels[h], " (", round(estGroupSize[h], 4), ")" ), xlab="", ylab="", text.size=1.3, 
        zlab=expression(p[ij]), label=FALSE, show.margins=FALSE,   label.lines=FALSE, sorted=FALSE, cum.margins=FALSE, colmar=0.5, rowmar=0.5, hide.duplicates=FALSE, 
        dotchar = 21, dotcolor = "steelblue", bg = "skyblue", dotsize = 1.4*sqrt(max(gewTransProb[,,h])*scal)*2/max(strwidth(19), strheight(19)) ) )
    
    if (plotPooled) {          
    
        dev.new(width=4,height=4) 
        par(mai=c(0.2, 0.2, 0.2, 0.25)) 
        balloonplot( k, j, c(t(relNjk)), 
            scale.method = c("volume"), main="pooled", xlab="", ylab="", text.size=1.3, 
            zlab=expression(p[ij]), label=FALSE, show.margins=FALSE, label.lines=FALSE, sorted=FALSE, cum.margins=FALSE, colmar=0.5, rowmar=0.5, hide.duplicates=FALSE, 
            dotchar = 21, dotcolor = "steelblue", bg ="skyblue", dotsize = 1.4*sqrt(max(relNjk)*scal) *2/max(strwidth(19), strheight(19))   )    
    }    
    
    if (plotContTable) {    
            
        Njk.h <- array(0,c(outList$K+1,outList$K+1,outList$Prior$H),dimnames=dimnames(outList$Njk.i) )
        for (h in 1:outList$Prior$H) Njk.h[,,h] <- apply(outList$Njk.i[,,class==h], c(1, 2), sum)
        tabRowSums <- as.matrix( apply(Njk.h, c(1,3), sum) )
        colnames(tabRowSums) <- grLabels               
        
        if (printContTable) {
            
            print( xtable( rbind(tabRowSums, sum=colSums(tabRowSums)), caption="Row sums of the contingency table (absolute transition frequencies) for each group.", digits=rep(0, 1+outList$Prior$H ) ) )
            cat("\n")
            print( xtable( round(100*rbind(tabRowSums, sum=colSums(tabRowSums))/matrix(colSums(tabRowSums), outList$K+2, outList$Prior$H, byrow=TRUE), 2), 
                    caption="'Relative' row sums (in percent) of the contingency table (absolute transition frequencies) for each group.", digits=rep(2, 1+outList$Prior$H ) ) ) 
            cat("\n")
        }
     
        relRowSums <- tabRowSums/matrix(colSums(tabRowSums), outList$K+1, outList$Prior$H, byrow=TRUE)     
        relTransFreq <- estTransProb * 0
        for (h in 1:outList$Prior$H) {    
            relTransFreq[,,h] <- estTransProb[,,h] * relRowSums[,h] 
        }
        Njk <- apply(outList$Njk.i, c(1, 2), sum)         
        relNjkMat <- Njk/sum(Njk)           
        scalW <- 1.0/max(relTransFreq, relNjkMat)        
        k   <- rep(as.numeric(dimnames( outList$Data$dataFile)[[1]]), outList$K+1)
        j   <- rep(as.numeric(dimnames( outList$Data$dataFile)[[1]]), each = outList$K+1)
                
        dev.new()
        par(mfrow=c(ceiling(outList$Prior$H/2),2), mai=c(0.2, 0.2, 0.2, 0.25))               
        sapply(1:outList$Prior$H, function(h) balloonplot( k, j, c(t(relTransFreq[,,h])), 
            scale.method = c("volume"), main=paste(grLabels[h], " (", round(estGroupSize[h], 4), ")" ), xlab="", ylab="", text.size=1.3, 
            zlab=expression(p[ij]), label=FALSE, show.margins=FALSE,   label.lines=FALSE, sorted=FALSE, cum.margins=FALSE, colmar=0.5, rowmar=0.5, hide.duplicates=FALSE, 
            dotchar = 21, dotcolor = "darkgreen", bg = "lightgreen", dotsize = 1.4*sqrt(max(relTransFreq[,,h])*scalW)*2/max(strwidth(19), strheight(19)) ))
            
        if (plotContPooled) {
        
            dev.new(width=4,height=4)
            par(mai=c(0.2, 0.2, 0.2, 0.25)) 
            balloonplot( k, j, c(t(relNjkMat)), 
                scale.method = c("volume"), main="pooled", xlab="", ylab="", text.size=1.3, 
                zlab=expression(p[ij]), label=FALSE, show.margins=FALSE, label.lines=FALSE, sorted=FALSE, cum.margins=FALSE, colmar=0.5, rowmar=0.5, hide.duplicates=FALSE, 
                dotchar = 21, dotcolor = "darkgreen", bg ="lightgreen", dotsize = 1.4*sqrt(max(relNjkMat)*scalW) *2/max(strwidth(19), strheight(19))   )    
        }    
    
    }
    
    transList <- list( relNjk=relNjk, contTable=tabRowSums, relTransFreq=relTransFreq, relNjkMat=relNjkMat )
    return( invisible( transList ) )
}
