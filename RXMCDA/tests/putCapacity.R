library(kappalab)
library(RXMCDA)

tree = newXMLDoc()

newXMLNode("xmcda:XMCDA", namespace = c("xsi" = "http://www.w3.org/2001/XMLSchema-instance", "xmcda" = "http://www.decision-deck.org/2009/XMCDA-2.1.0"), parent=tree)

mu<-capacity(0:1023)

a <- Mobius(mu)

critIDs <- c("g1","g2","g3","g4","g5","g6","g7","g8","g9","g10")

putCapacity(tree, a, critIDs, mcdaConcept="capacity")

capa<-getMobiusCapacities(tree, critIDs,10, 10, mcdaConcept="capacity")

stopifnot(identical(a,capa[[1]]))

