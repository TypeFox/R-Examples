library(MCDA)

# the separation threshold

epsilon <-0.05

# the performance table

performanceTable <- rbind(
  c(3,10,1),
  c(4,20,2),
  c(2,20,0),
  c(6,40,0),
  c(30,30,3))

rownames(performanceTable) <- c("RER","METRO1","METRO2","BUS","TAXI")

colnames(performanceTable) <- c("Price","Time","Comfort")

# ranks of the alternatives

alternativesAssignments <- c("good","medium","medium","bad","bad")

names(alternativesAssignments) <- row.names(performanceTable)

# criteria to minimize or maximize

criteriaMinMax <- c("min","min","max")

names(criteriaMinMax) <- colnames(performanceTable)

# number of break points for each criterion

criteriaNumberOfBreakPoints <- c(3,4,4)

names(criteriaNumberOfBreakPoints) <- colnames(performanceTable)

# ranks of the categories

categoriesRanks <- c(1,2,3)

names(categoriesRanks) <- c("good","medium","bad")

x<-UTADIS(performanceTable, criteriaMinMax, criteriaNumberOfBreakPoints, alternativesAssignments, categoriesRanks,0.1)

# filtering out category "good" and assigment examples "RER" and "TAXI" 

y<-UTADIS(performanceTable, criteriaMinMax, criteriaNumberOfBreakPoints, alternativesAssignments, categoriesRanks,0.1, categoriesIDs=c("medium","bad"), alternativesIDs=c("METRO1","METRO2","BUS"))

# working furthermore on only 2 criteria : "Comfort" and "Time"

z<-UTADIS(performanceTable, criteriaMinMax, criteriaNumberOfBreakPoints, alternativesAssignments, categoriesRanks,0.1, criteriaIDs=c("Comfort","Time"))

stopifnot(x$optimum ==0 && y$optimum ==0 && z$optimum ==0)



