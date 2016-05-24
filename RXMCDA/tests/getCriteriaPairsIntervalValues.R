library(RXMCDA)

tree <- xmlTreeParse(system.file("extdata","testFile.xml",package="RXMCDA"), useInternalNodes=TRUE)

critIDs <- getCriteriaIDs(tree)

pairsVals <- getCriteriaPairsIntervalValues(tree, critIDs[[1]],mcdaConcept="interactionIntervals")

stopifnot(all.equal(pairsVals[[1]], rbind(c(5,3,0.5,1))))
