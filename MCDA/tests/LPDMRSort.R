# ranking some students (from the article on Integrating Large Performance Differences in MR Sort by Meyer and Olteanu, 2015)

library(MCDA)

# the performance table

performanceTable <- rbind(c(10,10,9),c(10,9,10),c(9,10,10),c(9,9,10),c(9,10,9),c(10,9,9),
                          c(10,10,7),c(10,7,10),c(7,10,10),c(9,9,17),c(9,17,9),c(17,9,9),
                          c(7,10,17),c(10,17,7),c(17,7,10),c(7,17,10),c(17,10,7),c(10,7,17),
                          c(7,9,17),c(9,17,7),c(17,7,9),c(7,17,9),c(17,9,7),c(9,7,17))

profilesPerformances <- rbind(c(10,10,10),c(0,0,0))

vetoPerformances <- rbind(c(7,7,7),c(0,0,0))

dictatorPerformances <- rbind(c(17,17,17),c(0,0,0))

rownames(performanceTable) <- c("a1","a2","a3","a4","a5","a6","a7","a8","a9","a10","a11","a12","a13","a14","a15","a16","a17","a18","a19","a20","a21","a22","a23","a24")

rownames(profilesPerformances) <- c("P","F")

rownames(vetoPerformances) <- c("P","F")

rownames(dictatorPerformances) <- c("P","F")

colnames(performanceTable) <- c("c1","c2","c3")

colnames(profilesPerformances) <- c("c1","c2","c3")

colnames(vetoPerformances) <- c("c1","c2","c3")

colnames(dictatorPerformances) <- c("c1","c2","c3")

lambda <- 0.5

weights <- c(1/3,1/3,1/3)

names(weights) <- c("c1","c2","c3")

categoriesRanks <-c(1,2)

names(categoriesRanks) <- c("P","F")

criteriaMinMax <- c("max","max","max")

names(criteriaMinMax) <- colnames(performanceTable)

assignments <-rbind(c("P","P","P","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F"),
                    c("P","P","P","F","F","F","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P"),
                    c("P","P","P","F","F","F","F","F","F","F","F","F","P","P","P","P","P","P","F","F","F","F","F","F"),
                    c("P","P","P","F","F","F","P","P","P","P","P","P","P","P","P","P","P","P","F","F","F","F","F","F"),
                    c("P","P","P","F","F","F","F","F","F","P","P","P","F","F","F","F","F","F","F","F","F","F","F","F"),
                    c("P","P","P","F","F","F","F","F","F","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P"),
                    c("P","P","P","F","F","F","F","F","F","P","P","P","P","P","P","P","P","P","F","F","F","F","F","F"))

colnames(assignments) <- rownames(performanceTable)

majorityRules <- c("V","D","v","d","dV","Dv","dv")

for(i in 1:7)
{
  ElectreAssignments<-LPDMRSort(performanceTable, profilesPerformances, 
                                weights, criteriaMinMax, lambda, criteriaVetos=vetoPerformances, criteriaDictators=dictatorPerformances, majorityRule = majorityRules[i])

  print(ElectreAssignments)
  print(all(ElectreAssignments == assignments[i,]))

  stopifnot(all(assignments[i,]== ElectreAssignments))
}