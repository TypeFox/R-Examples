library(RXMCDA)

tree <- xmlTreeParse(system.file("extdata", "testFile.xml", package="RXMCDA"), useInternalNodes=TRUE)
categoriesIDs <- getCategoriesIDs(tree)
intervalValues <- getCategoriesIntervalValues(tree, categoriesIDs[[1]])

stopifnot(all.equal(intervalValues[[1]], matrix(data = c(1, 4, 2, 1, 6, NA), ncol = 3)))
