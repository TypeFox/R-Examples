library(RXMCDA)

alternativesIDs <-  c("a01", "a02", "a03", "a04")
categoriesIDs <- c("c01", "c02", "c03", "c04")
altAff = rbind(c(FALSE,TRUE,TRUE,TRUE), c(FALSE,TRUE,FALSE,FALSE), c(TRUE,TRUE,TRUE,TRUE), c(TRUE,TRUE,TRUE,FALSE))

tree = newXMLDoc()
newXMLNode("xmcda:XMCDA", namespace = c("xsi" = "http://www.w3.org/2001/XMLSchema-instance", "xmcda" = "http://www.decision-deck.org/2009/XMCDA-2.1.0"), parent=tree)

putAlternativesAffectations(tree, altAff, alternativesIDs, categoriesIDs, TRUE)

stopifnot(all.equal(getAlternativesAffectations(tree, alternativesIDs, categoriesIDs)[[1]], altAff))