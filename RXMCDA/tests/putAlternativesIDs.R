library(RXMCDA)

altIDs <- c("x","y","z")

tree = newXMLDoc()

newXMLNode("xmcda:XMCDA", namespace = c("xsi" = "http://www.w3.org/2001/XMLSchema-instance", "xmcda" = "http://www.decision-deck.org/2009/XMCDA-2.1.0"), parent=tree)

putAlternativesIDs(tree,altIDs)

altIDs2 <- getAlternativesIDs(tree)

stopifnot(identical(altIDs,altIDs2[[1]]))
