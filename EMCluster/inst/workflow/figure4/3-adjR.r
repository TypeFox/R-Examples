rm(list = ls())

library(EMCluster)
# source("./R/information.r")
load("./data/simu.data.rda")
# load("./data/simu.em.rda")
load("./data/simu.RndEM.rda")

for(k0 in 1:6){
  emobj0 <- ret.save[[k0]]

  cat("K = ", k0, "\n", sep = "")
  print(RRand(emobj0$class, as.integer(da$id)))
}
