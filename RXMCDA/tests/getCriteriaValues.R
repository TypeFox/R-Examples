library(RXMCDA)

tree <- xmlTreeParse(system.file("extdata","testFile.xml",package="RXMCDA"), useInternalNodes=TRUE)

critIDs <- getCriteriaIDs(tree)

critVals <- getCriteriaValues(tree, critIDs[[1]])

stopifnot(all.equal(critVals[[1]], rbind(c(5,0.11), c(3,0.173363501599438), c(2,0.439744089880627))))
