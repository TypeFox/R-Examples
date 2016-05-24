library(RXMCDA)

tree <- xmlTreeParse(system.file("extdata","testFile.xml",package="RXMCDA"), useInternalNodes=TRUE)

tables <- getPerformanceTables(tree)

stopifnot(all.equal(dim(tables[[1]]), c(15,5)))


