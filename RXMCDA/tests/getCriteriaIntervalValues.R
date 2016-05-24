library(RXMCDA)

tree <- xmlTreeParse(system.file("extdata","testFile.xml",package="RXMCDA"), useInternalNodes=TRUE)

critIDs <- getCriteriaIDs(tree)

intVals <- getCriteriaIntervalValues(tree, critIDs[[1]])

stopifnot(all.equal(intVals[[1]], t(as.matrix(c(5,0.5,1)))))
