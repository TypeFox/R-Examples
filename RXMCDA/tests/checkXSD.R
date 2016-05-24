library(RXMCDA)

tree = newXMLDoc()

newXMLNode("xmcda:XMCDA", namespace = c("xsi" = "http://www.w3.org/2001/XMLSchema-instance", "xmcda" = "http://www.decision-deck.org/2009/XMCDA-2.0.0"), parent=tree)

root<-xmlRoot(tree)

criteria<-newXMLNode("criteria", parent=root, namespace=c())

criterion<-newXMLNode("criterion",attrs = c(id="g1"), parent=criteria, namespace=c())

y<-checkXSD(tree)

stopifnot(y[1] == 1)

tree = newXMLDoc()

newXMLNode("xmcda:XMCDA", namespace = c("xsi" = "http://www.w3.org/2001/XMLSchema-instance", "xmcda" = "http://www.decision-deck.org/2009/XMCDA-2.1.0"), parent=tree)

root<-xmlRoot(tree)

criteria<-newXMLNode("criteria", parent=root, namespace=c())

criterion<-newXMLNode("criterion",attrs = c(id="g1"), parent=criteria, namespace=c())

y<-checkXSD(tree)

stopifnot(y[1] == 1)

tree = newXMLDoc()

newXMLNode("xmcda:XMCDA", namespace = c("xsi" = "http://www.w3.org/2001/XMLSchema-instance", "xmcda" = "http://www.decision-deck.org/2012/XMCDA-2.2.0"), parent=tree)

root<-xmlRoot(tree)

criteria<-newXMLNode("criteria", parent=root, namespace=c())

criterion<-newXMLNode("criterion",attrs = c(id="g1"), parent=criteria, namespace=c())

y<-checkXSD(tree)

stopifnot(y[1] == 1)

