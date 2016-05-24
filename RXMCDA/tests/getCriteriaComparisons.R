library(RXMCDA)

tree <- xmlTreeParse(system.file("extdata","testFile.xml",package="RXMCDA"), useInternalNodes=TRUE)

critIDs <- getCriteriaIDs(tree)

comps <- getCriteriaComparisons(tree, critIDs[[1]])

stopifnot(all.equal(dim(comps[[1]]), c(2,3)))


