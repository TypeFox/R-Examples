library(MCDA)

performanceTable <- matrix(runif(6*9), ncol=6)

row.names(performanceTable) <- c("x1","x2","x3","x4","x5","x6","x7","x8","x9")

colnames(performanceTable) <- c("g1","g2","g3","g4","g5","g6")

criteriaMinMax <- c("min","max","min","max","min","max")

names(criteriaMinMax) <- c("g1","g2","g3","g4","g5","g6")

# plotRadarPerformanceTable(performanceTable, criteriaMinMax, overlay=TRUE)

plotRadarPerformanceTable(performanceTable, criteriaMinMax, alternativesIDs <- c("x1","x2","x3","x4"), criteriaIDs<-c("g1","g3","g4","g5","g6"), overlay=FALSE)

# plotRadarPerformanceTable(performanceTable, criteriaMinMax, alternativesIDs <- c("x1","x2"), criteriaIDs<-c("g1","g3","g4","g5","g6"), overlay=FALSE)

