library(RXMCDA)

tree = newXMLDoc()

newXMLNode("xmcda:XMCDA", namespace = c("xsi" = "http://www.w3.org/2001/XMLSchema-instance", "xmcda" = "http://www.decision-deck.org/2009/XMCDA-2.1.0"), parent=tree)

critIDs <- c("g1","g2","g3","g4")

pairsVals <- rbind(c(1,2,0.17),c(2,3,0.5), c(3,4,0.16))

putCriteriaPairsValues(tree,pairsVals,critIDs)

pairsVals2 <- getCriteriaPairsValues(tree,critIDs)

stopifnot(identical(pairsVals,pairsVals2[[1]]))
