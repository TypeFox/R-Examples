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

# lower profiles of the categories (best category in the first position of the list)

categoriesLowerProfiles <- rbind(c(3, 11, 3),c(7, 25, 2),c(NA,NA,NA))

colnames(categoriesLowerProfiles) <- colnames(performanceTable)

rownames(categoriesLowerProfiles)<-c("Good","Medium","Bad")

# criteria to minimize or maximize

criteriaMinMax <- c("min","min","max")

names(criteriaMinMax) <- colnames(performanceTable)

# vetos

criteriaVetos <- rbind(c(10, NA, NA),c(NA, NA, 1),c(NA,NA,NA))

colnames(criteriaVetos) <- colnames(performanceTable)
rownames(criteriaVetos) <- c("Good","Medium","Bad")

# weights

criteriaWeights <- c(1,3,2)

names(criteriaWeights) <- colnames(performanceTable)


# MRSort

assignments<-MRSort(performanceTable, categoriesLowerProfiles, criteriaWeights, criteriaMinMax, 3, criteriaVetos = criteriaVetos)

stopifnot(all(assignments == c("Good","Medium","Bad","Bad","Bad")))

# un peu de filtrage

assignments<-MRSort(performanceTable, categoriesLowerProfiles, criteriaWeights, criteriaMinMax, 2, categoriesIDs = c("Medium","Bad"), criteriaIDs = c("Price","Time"), alternativesIDs = c("RER", "BUS"))

stopifnot(all(assignments == c("Medium","Bad")))

# un test pour combiner tous les filtrages avec le veto

assignments<-MRSort(performanceTable, categoriesLowerProfiles, criteriaWeights, criteriaMinMax, 2, criteriaVetos = criteriaVetos, categoriesIDs = c("Medium","Bad"), criteriaIDs = c("Price","Time"), alternativesIDs = c("RER", "BUS"))

stopifnot(all(assignments == c("Medium","Bad")))
