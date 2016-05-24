rm(list = ls())

library(EMCluster)
# source("./R/information.r")
# load("./data/iris.em.rda")
load("./data/iris.RndEM.rda")

id <- as.integer(iris[, 5])

for(k0 in 1:7){
  emobj0 <- ret.save[[k0]]

  cat("K = ", k0, "\n", sep = "")
  print(RRand(emobj0$class, id))
}
