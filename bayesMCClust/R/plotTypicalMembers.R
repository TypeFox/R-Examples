plotTypicalMembers <-
function(outList, myObsList, classProbs, noTypMemb=7, moreTypMemb=c(10,25,50,100,200,500,1000), grLabels=paste("Group", 1:outList$Prior$H)) {

    stopifnot( noTypMemb==length(moreTypMemb) )

    stopifnot( length(myObsList)==outList$N )
    stopifnot( nrow(classProbs)==outList$N & ncol(classProbs)==outList$Prior$H ) 

    typicalMemb  <- matrix(0, outList$Prior$H, noTypMemb)
    typicalMemb2 <- matrix(0, outList$Prior$H, noTypMemb)

    for (h in 1:outList$Prior$H) {  
        typMembTemp <- sort.int( classProbs[,h], decreasing = TRUE, index.return=TRUE)
        typicalMemb[h,]  <- typMembTemp$ix[1:noTypMemb]  
        typicalMemb2[h,] <- typMembTemp$ix[moreTypMemb] 
    }

    dev.new(width=12,height=8) 
    par(mfrow=c(outList$Prior$H, noTypMemb), mai=c(0.6, 0.3, 0.35, 0.1) ) # c(bottom, left, top, right)
    
    for (h in 1:outList$Prior$H) {
        for (i in 1:noTypMemb) { 
            plot(myObsList[[typicalMemb[h,i]]], type="p", ylab="k", bty="l", xlab=if (h==outList$Prior$H) paste("most typical group memb no", i, sep=" ") else " ", 
                    ylim=range(as.numeric(dimnames( outList$Data$dataFile)[[1]])), pch=19, bg="black", las=1, yaxt="n", 
                    cex=1.3, cex.axis=1.7, cex.main = 1.9, main=if ((i+2)%%6==0) grLabels[h] )  
            axis(2, at = as.numeric(dimnames( outList$Data$dataFile)[[1]]), labels = dimnames(outList$Data$dataFile)[[1]], cex.axis=1.7, las=1) # c("K","0","1","2","3")
            
        }
    }

    dev.new(width=12,height=8) 
    par(mfrow=c(outList$Prior$H, noTypMemb), mai=c(0.6, 0.3, 0.35, 0.1) ) # c(bottom, left, top, right)

    for (h in 1:outList$Prior$H) {
        for (i in 1:noTypMemb) { 
            plot(myObsList[[typicalMemb2[h,i]]], type="p", ylab="k", bty="l", xlab=if (h==outList$Prior$H) paste("most typ gr memb no", moreTypMemb[i], sep=" ") else " ", 
                    ylim=range(as.numeric(dimnames( outList$Data$dataFile)[[1]])), pch=19, bg="black", las=1, yaxt="n", # 
                    cex=1.3, cex.axis=1.7, cex.main = 1.9, main=if ((i+2)%%6==0) grLabels[h] )  
            axis(2, at = as.numeric(dimnames( outList$Data$dataFile)[[1]]), labels = dimnames(outList$Data$dataFile)[[1]], cex.axis=1.7, las=1) # c("K","0","1","2","3")
        }
    }
    
    return( invisible( list( typicalMemb=typicalMemb, typicalMemb2=typicalMemb2 ) ) )

}
