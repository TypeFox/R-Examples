library(RXMCDA)

tree <- xmlTreeParse(system.file("extdata","testFile.xml",package="RXMCDA"), useInternalNodes=TRUE)

critIDs <- getCriteriaIDs(tree)

pairsVals <- getCriteriaPairsValues(tree, critIDs[[1]],mcdaConcept="interactionValues")

stopifnot(all.equal(pairsVals[[1]], rbind(c(5,3,0.17),c(5,2,-0.238940469688695), c(5,4,0.160372942186159))))
