library(MCDA)

## This test example is the same as http://www.tandfonline.com/doi/abs/10.1080/09537287.2013.798706

performanceTableMin <- t(matrix(c(78,87,79,19,8,68,74,8,90,89,74.5,9,20,81,30),nrow=3,ncol=5, byrow=TRUE)) 
performanceTable <- t(matrix(c(80,87,86,19,8,70,74,10,90,89,75,9,33,82,30),nrow=3,ncol=5, byrow=TRUE))
performanceTableMax <- t(matrix(c(81,87,95,19,8,72,74,15,90,89,75.5,9,36,84,30),nrow=3,ncol=5, byrow=TRUE))  

row.names(performanceTable) <- c("Yield","Toxicity","Cost","Separation","Odour")
colnames(performanceTable) <- c("Route One","Route Two","Route Three")
row.names(performanceTableMin) <- row.names(performanceTable)
colnames(performanceTableMin) <- colnames(performanceTable)
row.names(performanceTableMax) <- row.names(performanceTable)
colnames(performanceTableMax) <- colnames(performanceTable)

weights <- c(0.339,0.077,0.434,0.127,0.023) 
names(weights) <- row.names(performanceTable)

criteriaMinMax <- c("max", "max", "max", "max", "max")
names(criteriaMinMax) <- row.names(performanceTable)

overall1 <- MARE(performanceTableMin, performanceTable, performanceTableMax, weights, criteriaMinMax)

overall2 <- MARE(performanceTableMin, performanceTable, performanceTableMax, weights, criteriaMinMax, alternativesIDs = alternativesIDs <- c("Route Two","Route Three"), criteriaIDs = criteriaIDs <- c("Yield","Toxicity","Cost","Separation"))

s1 = structure(c(0.7932, 0.8336, 0.8789, 0.5366, 0.5541, 0.5854, 0.5332, 0.5961, 0.6147), .Dim = c(3L, 3L), .Dimnames = list(c("Minimum", "Most Likely", "Maximum"), c("Route One", "Route Two", "Route Three")))

stopifnot(round(overall1,4) == s1)

s2 = structure(c(0.6058, 0.6389, 0.7081, 0.6993, 0.8597, 0.9009), .Dim = c(3L, 2L), .Dimnames = list(c("Minimum", "Most Likely", "Maximum"), c("Route Two", "Route Three")))

stopifnot(round(overall2,4) == s2)
