library(MCDA)

performanceTable <- matrix(runif(5*9), ncol=5)

row.names(performanceTable) <- c("x1","x2","x3","x4","x5","x6","x7","x8","x9")

colnames(performanceTable) <- c("g1","g2","g3","g4", "g5")

normalizationTypes <- c("percentageOfMax","rescaling","standardization","scaleToUnitLength", "none")

names(normalizationTypes) <- c("g1","g2","g3","g4","g5")

normalizedPerformanceTable <- normalizePerformanceTable(performanceTable,normalizationTypes)

stopifnot(
  (max(normalizedPerformanceTable[,"g1"]) == 1) 
  & (max(normalizedPerformanceTable[,"g2"])-min(normalizedPerformanceTable[,"g2"])==1) 
  & (abs(sd(normalizedPerformanceTable[,"g3"])-1) < 0.001) 
  & (abs(sqrt(sum(normalizedPerformanceTable[,"g4"]^2)) - 1) < 0.001)
  )

normalizedPerformanceTable <- normalizePerformanceTable(performanceTable,normalizationTypes,criteriaIDs=c("g2","g3"))

stopifnot(
  (max(normalizedPerformanceTable[,"g2"])-min(normalizedPerformanceTable[,"g2"])==1) 
  & (abs(sd(normalizedPerformanceTable[,"g3"])-1) < 0.001)
  )

normalizedPerformanceTable <- normalizePerformanceTable(performanceTable,normalizationTypes,alternativesIDs=c("x2","x3"))

stopifnot(
  (max(normalizedPerformanceTable[,"g1"]) == 1) 
  & (max(normalizedPerformanceTable[,"g2"])-min(normalizedPerformanceTable[,"g2"])==1) 
  & (abs(sd(normalizedPerformanceTable[,"g3"])-1) < 0.001) 
  & (abs(sqrt(sum(normalizedPerformanceTable[,"g4"]^2)) - 1) < 0.001)
)

normalizedPerformanceTable <- normalizePerformanceTable(performanceTable,normalizationTypes,alternativesIDs=c("x2","x3"),criteriaIDs=c("g1","g2"))

stopifnot(
  (max(normalizedPerformanceTable[,"g1"]) == 1) 
  & (max(normalizedPerformanceTable[,"g2"])-min(normalizedPerformanceTable[,"g2"])==1) 
)