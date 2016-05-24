library(RXMCDA)

tree <- xmlTreeParse(system.file("extdata", "testFile.xml", package="RXMCDA"), useInternalNodes=TRUE)
categoriesIDs <- getCategoriesIDs(tree)
categoriesValues <- getCategoriesValues(tree, categoriesIDs[[1]])

stopifnot(all.equal(categoriesValues[[1]], matrix(data = c(1, 2, 100, -100), ncol = 2)))