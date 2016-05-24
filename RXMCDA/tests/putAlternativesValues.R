library(RXMCDA)

altIDs <- c("x","y","z")

altVals <- rbind(c(1,1),c(2,0.5),c(3,0.2))

tree = newXMLDoc()

newXMLNode("xmcda:XMCDA", namespace = c("xsi" = "http://www.w3.org/2001/XMLSchema-instance", "xmcda" = "http://www.decision-deck.org/2009/XMCDA-2.1.0"), parent=tree)

putAlternativesValues(tree,altVals,altIDs)

altVals2 <- getAlternativesValues(tree, altIDs)

stopifnot(identical(altVals,altVals2[[1]]))
