library(RXMCDA)

tree = newXMLDoc()

newXMLNode("xmcda:XMCDA", 
           namespace = c("xsi" = "http://www.w3.org/2001/XMLSchema-instance", 
           "xmcda" = "http://www.decision-deck.org/2009/XMCDA-2.1.0"), 
           parent=tree)

root <- getNodeSet(tree, "/xmcda:XMCDA")
altComp <- newXMLNode("alternativesComparisons", parent=root[[1]] , namespace=c())
pairs <- newXMLNode("pairs", parent=altComp, namespace=c())


pair <- newXMLNode("pair", parent=pairs, namespace=c())
initial <- newXMLNode("initial", parent=pair)
newXMLNode("alternativeID", "a01", parent=initial, namespace=c())
terminal <- newXMLNode("terminal", parent=pair, namespace=c())
newXMLNode("alternativeID", "a02", parent=terminal, namespace=c())
value <- newXMLNode("value", parent=pair, namespace=c())
newXMLNode("real", "1", parent=value, namespace=c())

pair <- newXMLNode("pair", parent=pairs, namespace=c())
initial <- newXMLNode("initial", parent=pair)
newXMLNode("alternativeID", "a01", parent=initial, namespace=c())
terminal <- newXMLNode("terminal", parent=pair, namespace=c())
newXMLNode("alternativeID", "a03", parent=terminal, namespace=c())
value <- newXMLNode("value", parent=pair, namespace=c())
newXMLNode("real", "9", parent=value, namespace=c())

alternativesIDs <- c("a01", "a02", "a03")

altCompVal <- getAlternativesComparisonsValues(tree, alternativesIDs)

stopifnot(nrow(altCompVal[[1]]) == 2)
stopifnot(all.equal(altCompVal[[1]][1, ], c(1, 2, 1)))
stopifnot(all.equal(altCompVal[[1]][2, ], c(1, 3, 9)))