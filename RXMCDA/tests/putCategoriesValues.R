library(RXMCDA)

categoriesIDs <- c("c01", "c02", "c03", "c04")
categoriesValues <- rbind(c(1, 0.4), c(2, 0.5), c(4, 0.2))

tree = newXMLDoc()
newXMLNode("xmcda:XMCDA", namespace = c("xsi" = "http://www.w3.org/2001/XMLSchema-instance", "xmcda" = "http://www.decision-deck.org/2009/XMCDA-2.1.0"), parent=tree)

putCategoriesValues(tree, categoriesValues, categoriesIDs)

stopifnot(all.equal(getCategoriesValues(tree, categoriesIDs)[[1]], categoriesValues))