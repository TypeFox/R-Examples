geneNetworkSummary <-
function(ARTIVAnet, edgesThreshold){

 # Only interactions with a posterior probability higher than the specified threshold are written
 ARTIVAsubRes = ARTIVAnet[ARTIVAnet$PostProb >= edgesThreshold,]

 # Interactions are ordered 
 ARTIVAsubRes = ARTIVAsubRes[order(ARTIVAsubRes$PostProb, decreasing = T),]
    
 # To get the signs of the interactions
 InteractionSigns = rep("NA", nrow(ARTIVAsubRes))
 InteractionSigns[ARTIVAsubRes$CoeffMean >= 0] = "+"
 InteractionSigns[ARTIVAsubRes$CoeffMean < 0]  = "-"

 # Information to be written
 ResTable =  data.frame(cbind(ARTIVAsubRes$Parent, ARTIVAsubRes$Target, ARTIVAsubRes$PostProb,
                                 ARTIVAsubRes$CPstart, ARTIVAsubRes$CPend, InteractionSigns))  
 colnames(ResTable) <- c("parentGene", "targetGene", "postProb", "CPstart", 
                            "CPend", "interactionSign")
 #Calculate the number of pages (20 lines per page)
lineNumber = 20
pageNumber = floor(nrow(ResTable)/lineNumber) + 1

j = 1
for(i in 1:pageNumber){
    textplot(ResTable[j:(j + lineNumber - 1),], cex = 1, cmar = 1.1, rmar = 1, 
             show.rownames = F, show.colnames = T,
             halign = "center", valign = "top")
    title("ARTIVA summary page\n(interactions are arranged in order of decreasing confidence level)")
    
    j = j + lineNumber
# end of for()
} 

# last page
if(nrow(ResTable) > (pageNumber * lineNumber)){

    textplot(ResTable[(pageNumber * lineNumber) + 1:nrow(ResTable),], 
             cex = 1, cmar = 1.1, rmar = 1, 
             show.rownames = F, show.colnames = T,
             halign = "center", valign = "top")
    title("ARTIVA summary page\n(interactions are arranged in order of decreasing confidence level)")
    
# end of if()
}

# End of function
}
