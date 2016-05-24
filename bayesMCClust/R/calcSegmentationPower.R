calcSegmentationPower <-
function(outList, classProbs, class, printXtable=TRUE, calcSharp=TRUE, printSharpXtable=TRUE, grLabels=paste("Group", 1:outList$Prior$H) ) {

    maxProbs <- apply(classProbs, 1, max)  

    segPowTab <- matrix(0,outList$Prior$H+1, 6, dimnames=list(c(grLabels,"overall"), c(  names( summary(maxProbs))) ))

    for (h in 1:outList$Prior$H) segPowTab[h,] <- c( summary(maxProbs[class==h]))

    segPowTab[outList$Prior$H+1,] <- c( summary(maxProbs))

    if ( printXtable ) {
    
        print( xtable( segPowTab, digits=4, caption=paste("Segmentation power: reported are some summary statistics for the maximum individual posterior classification probabilities",
                                                          "for all individuals within a certain cluster as well as for all individuals."), label="tab:segPower" ) )        
        cat("\n")     
    }
    
    if ( calcSharp ) {
    
        sharp <- apply(classProbs, 1, function(x) max(x) - sort(x, decreasing = TRUE)[2]) 
    
        sharpTab <- matrix(0,outList$Prior$H+1, 6, dimnames=list(c(grLabels,"overall"), c(  names( summary(sharp))) ))
    
        for (h in 1:outList$Prior$H) sharpTab[h,] <- c( summary(sharp[class==h]))
        
        sharpTab[outList$Prior$H+1,] <- c( summary(sharp))
    
        if ( printXtable ) {
    
            print( xtable( sharpTab, digits=4, caption="'Sharpness': reported are some summary statistics for the difference between highest and second highest 
                                individual posterior classification probabilities within groups and overall", 
                           label="tab:sharpness" ) )
            cat("\n")     
        }
    
    }    
    
    return( invisible( list(segPowTab=segPowTab, sharpTab=sharpTab, maxProbs=maxProbs, sharp=sharp) ) )

}
