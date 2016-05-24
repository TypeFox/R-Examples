calcEntropy <-
function(outList, classProbs, class, grLabels=paste("Group", 1:outList$Prior$H), printXtable=TRUE ) {
 
    entr <- numeric(outList$Prior$H+1)
    for (h in 1:outList$Prior$H) {
        entr[h] <- -sum(classProbs[,h]*log(classProbs[,h]))
    }
    entr[outList$Prior$H+1] <- -sum(classProbs*log(classProbs))

    entrTable <- cbind(entr, c(table(class), outList$N), entr/c(table(class), outList$N))

    dimnames(entrTable) <- list(c(grLabels, "total"), c("abs entr", "abs gr size", "rel entr"))
    
    if ( printXtable ) {
    
        print( xtable( entrTable, digits = c(0,2,0,2), caption="Contribution of each group to the total entropy (absolute and relative to group size -- based on ind post class probs)", 
                       label="tab:entropy"))
        
        cat("\n") 
    
    }
    
    return( invisible( entrTable ) )

}
