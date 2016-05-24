library(RXMCDA)

categoriesIDs <- c("c01", "c02", "c03", "c04")
categoriesIntervalValues <- rbind(c(1, 0.4, 0.7), c(2, 0.5, 0.5), c(4, 0.2, 0.9))

tree = newXMLDoc()
newXMLNode("xmcda:XMCDA", namespace = c("xsi" = "http://www.w3.org/2001/XMLSchema-instance", "xmcda" = "http://www.decision-deck.org/2009/XMCDA-2.1.0"), parent=tree)

putCategoriesIntervalValues(tree, categoriesIntervalValues, categoriesIDs)

stopifnot(all.equal(getCategoriesIntervalValues(tree, categoriesIDs)[[1]], categoriesIntervalValues))