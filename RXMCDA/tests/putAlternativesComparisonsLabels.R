library(RXMCDA)

tree = newXMLDoc()

newXMLNode("croquette:XMCDA", namespace = c("xsi" = "http://www.w3.org/2001/XMLSchema-instance", "croquette" = "http://www.decision-deck.org/2009/XMCDA-2.1.0"), parent=tree)

comps <- rbind(c("x", "y", "0.07"), c("y", "z", "0.01"))

altIDs <- c("x","y","z")

putAlternativesComparisonsLabels(tree,comps, mcdaConcept="newComparisons")

comps2 <- getAlternativesComparisonsLabels(tree, altIDs, mcdaConcept="newComparisons")

stopifnot(identical(comps,comps2[[1]]))
