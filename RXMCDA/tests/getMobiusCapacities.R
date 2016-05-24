library(RXMCDA)

tree <- xmlTreeParse(system.file("extdata","testFile.xml",package="RXMCDA"), useInternalNodes=TRUE)

critIDs <- getCriteriaIDs(tree)

capa <- getMobiusCapacities(tree, critIDs[[1]], 5, 5, mcdaConcept="mobiusCapacity")

stopifnot(class(capa[[1]])[[1]] == "Mobius.capacity")
