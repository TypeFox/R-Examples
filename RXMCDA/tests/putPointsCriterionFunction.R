library(RXMCDA)

tree = newXMLDoc()

newXMLNode("xmcda:XMCDA", namespace = c("xsi" = "http://www.w3.org/2001/XMLSchema-instance", "xmcda" = "http://www.decision-deck.org/2009/XMCDA-2.1.0"), parent=tree)

x<-list()
x<-c(x,list(g1=rbind(c(1,2),c(3,4))))
x<-c(x,list(g2=rbind(c(5,6),c(7,8),c(9,10))))
x<-c(x,list(g3=rbind(c(11,12))))
x<-c(x,list(g4=rbind(c(13,14),c(15,16))))

putPointsCriterionFunction(tree,x)

