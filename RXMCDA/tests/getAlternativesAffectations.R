library(RXMCDA)

tree <- xmlTreeParse(system.file("extdata", "testFile.xml", package="RXMCDA"), useInternalNodes=TRUE)

alternativesIDs <- getAlternativesIDs(tree)
categoriesIDs <- getCategoriesIDs(tree)
altAff <- getAlternativesAffectations(tree, alternativesIDs[[1]], categoriesIDs[[1]])

stopifnot(all.equal(altAff[[1]][1, ], c(FALSE, TRUE, TRUE, TRUE)))
stopifnot(all.equal(altAff[[1]][2, ], c(FALSE, TRUE, TRUE, FALSE)))
stopifnot(all.equal(altAff[[1]][3, ], c(FALSE, FALSE, TRUE, FALSE)))