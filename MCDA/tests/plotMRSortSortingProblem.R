library(MCDA)

# the performance table

performanceTable <- rbind(
  c(1,10,1),
  c(4,20,2),
  c(2,20,0),
  c(6,40,0),
  c(30,30,3))

rownames(performanceTable) <- c("RER","METRO1","METRO2","BUS","TAXI")

colnames(performanceTable) <- c("Price","Time","Comfort")

# lower profiles of the categories 
# (best category in the first position of the list)

categoriesLowerProfiles <- rbind(c(3, 11, 3),c(7, 25, 2),c(30,30,0))

colnames(categoriesLowerProfiles) <- colnames(performanceTable)

rownames(categoriesLowerProfiles)<-c("Good","Medium","Bad")

# criteria to minimize or maximize

criteriaMinMax <- c("min","min","max")

names(criteriaMinMax) <- colnames(performanceTable)

# lower bounds of the criteria to be used for the determination of value functions

criteriaLBs=c(0,5,0)

names(criteriaLBs) <- colnames(performanceTable)

# upper bounds of the criteria to be used for the determination of value functions

criteriaUBs=c(50,50,4)

names(criteriaUBs) <- colnames(performanceTable)

# weights

criteriaWeights <- c(1,3,2)

names(criteriaWeights) <- colnames(performanceTable)

assignments <- assignments<-MRSort(performanceTable, categoriesLowerProfiles, 
                                         criteriaWeights, criteriaMinMax, 3)

names(assignments) <- rownames(performanceTable)

plotMRSortSortingProblem(performanceTable, categoriesLowerProfiles, assignments, criteriaMinMax, criteriaUBs, criteriaLBs)