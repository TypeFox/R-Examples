library(RXMCDA)

tree = newXMLDoc()

newXMLNode("xmcda:XMCDA", 
           namespace = c("xsi" = "http://www.w3.org/2001/XMLSchema-instance", 
           "xmcda" = "http://www.decision-deck.org/2009/XMCDA-2.0.0"), 
           parent=tree)

root<-getNodeSet(tree, "/xmcda:XMCDA")

categories<-newXMLNode("categories", parent=root[[1]], namespace=c())

newXMLNode("category", attrs = c(id="c1"), parent=categories, namespace=c())
newXMLNode("category", attrs = c(id="c2"), parent=categories, namespace=c())

y<-getNodeSet(tree,"//categories")

stopifnot(getNumberOfCategories(y[[1]])[[1]] == 2)
