library(MCDA)

# the value functions

v<-list(
  Price = array(c(30, 0, 16, 0, 2, 0.0875), 
                dim=c(2,3), dimnames = list(c("x", "y"), NULL)), 
  Time = array(c(40, 0, 30, 0, 20, 0.025, 10, 0.9), 
               dim = c(2, 4), dimnames = list(c("x", "y"), NULL)), 
  Comfort = array(c(0, 0, 1, 0, 2, 0.0125, 3, 0.0125), 
                  dim = c(2, 4), dimnames = list(c("x", "y"), NULL)))

# the performance table

performanceTable <- rbind(
  c(3,10,1),
  c(4,20,2),
  c(2,20,0),
  c(6,40,0),
  c(30,30,3))

rownames(performanceTable) <- c("RER","METRO1","METRO2","BUS","TAXI")

colnames(performanceTable) <- c("Price","Time","Comfort")

# criteria to be minimized or maximized

criteriaMinMax <- c("min","min","max")

names(criteriaMinMax) <- colnames(performanceTable)

# the transformed performance table

normalizedPerformanceTable <- applyPiecewiseLinearValueFunctionsOnPerformanceTable(v,performanceTable)

correctResult <- rbind(c(0.08125,0.9,0),
                       c(0.075,0.025,0.0125),
                       c(0.0875,0.025,0),
                       c(0.0625,0,0),
                       c(0,0,0.0125))

rownames(correctResult) <- c("RER","METRO1","METRO2","BUS","TAXI")

colnames(correctResult) <- c("Price","Time","Comfort")

stopifnot(identical(normalizedPerformanceTable,correctResult))
