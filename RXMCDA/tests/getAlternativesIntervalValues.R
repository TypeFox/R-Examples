library(RXMCDA)

tree <- xmlTreeParse(system.file("extdata","testFile.xml",package="RXMCDA"), useInternalNodes=TRUE)

altIDs <- getAlternativesIDs(tree)

altVals <- getAlternativesIntervalValues(tree, altIDs[[1]], mcdaConcept="alternativesIntervalValues")

stopifnot(all.equal(dim(altVals[[1]]), c(2,3)))
