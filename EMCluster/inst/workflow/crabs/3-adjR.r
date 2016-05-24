rm(list = ls())

library(EMCluster)
library(MASS)
# source("./R/information.r")
# load("./data/crabs.em.rda")
load("./data/crabs.RndEM.rda")

crabs.id <- rep(1L, nrow(crabs))
crabs.id[crabs$sp == "B" & crabs$sex == "F"] <- 2L
crabs.id[crabs$sp == "O" & crabs$sex == "M"] <- 3L
crabs.id[crabs$sp == "O" & crabs$sex == "F"] <- 4L

for(k0 in 1:7){
  emobj0 <- ret.save[[k0]]

  cat("K = ", k0, "\n", sep = "")
  print(RRand(emobj0$class, crabs.id))
}
