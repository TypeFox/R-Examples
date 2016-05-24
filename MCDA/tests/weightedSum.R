library(MCDA)

performanceTable <- matrix(runif(3*4), ncol=3)

row.names(performanceTable) <- c("x1","x2","x3","x4")

colnames(performanceTable) <- c("g1","g2","g3")

weights <- c(1,2,3)

names(weights) <- c("g1","g2","g3")

overall1 <- weightedSum(performanceTable, weights)

overall2 <- weightedSum(performanceTable, weights, alternativesIDs <- c("x2","x3"), criteriaIDs <- c("g2","g3"))

stopifnot(length(overall2) == 2)


