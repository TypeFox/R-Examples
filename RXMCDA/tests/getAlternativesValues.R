library(RXMCDA)

tree <- xmlTreeParse(system.file("extdata","testFile.xml",package="RXMCDA"), useInternalNodes=TRUE)

altIDs <- getAlternativesIDs(tree)

altVals <- getAlternativesValues(tree, altIDs[[1]])

stopifnot(all.equal(dim(altVals[[1]]), c(15,2)))
