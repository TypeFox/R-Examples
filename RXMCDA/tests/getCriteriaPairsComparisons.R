library(RXMCDA)

tree <- xmlTreeParse(system.file("extdata","testFile.xml",package="RXMCDA"), useInternalNodes=TRUE)

critIDs <- getCriteriaIDs(tree)

comps <- getCriteriaPairsComparisons(tree, critIDs[[1]])

stopifnot(all.equal(dim(comps[[1]]), c(2,5)))

stopifnot(all.equal(comps[[1]], rbind(c(1,2,1,5,0.08),c(4,5,1,5,3))))
