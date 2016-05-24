library(RXMCDA)

critIDs <- c("x","y","z")

critVals <- rbind(c(1,1),c(2,0.5),c(3,0.2))

tree = newXMLDoc()

newXMLNode("xmcda:XMCDA", namespace = c("xsi" = "http://www.w3.org/2001/XMLSchema-instance", "xmcda" = "http://www.decision-deck.org/2009/XMCDA-2.1.0"), parent=tree)

putCriteriaValues(tree,critVals,critIDs)

critVals2 <- getCriteriaValues(tree, critIDs)

stopifnot(identical(critVals,critVals2[[1]]))
