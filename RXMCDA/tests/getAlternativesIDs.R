library(RXMCDA)

tree = newXMLDoc()

newXMLNode("xmcda:XMCDA", namespace = c("xsi" = "http://www.w3.org/2001/XMLSchema-instance", "xmcda" = "http://www.decision-deck.org/2009/XMCDA-2.0.0"), parent=tree)

root<-getNodeSet(tree, "/xmcda:XMCDA")

alternatives<-newXMLNode("alternatives", attrs=c(mcdaConcept="actions"), parent=root[[1]], namespace=c())

alternative<-newXMLNode("alternative",attrs = c(id="x1"), parent=alternatives, namespace=c())
alternative<-newXMLNode("alternative",attrs = c(id="x2"), parent=alternatives, namespace=c())
alternative<-newXMLNode("alternative",attrs = c(id="x3"), parent=alternatives, namespace=c())

y<-getNodeSet(tree,"//alternatives")

stopifnot(getAlternativesIDs(y[[1]])[[1]] == c("x1","x2","x3"))
