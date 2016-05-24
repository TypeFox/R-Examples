rm(list = ls())

library(MixSim)
set.seed(1234)

n <- 200
p <- 2
K <- 4
omega.bar <- 0.05
omega.check <- 0.10

param <- MixSim(BarOmega = omega.bar, MaxOmega = omega.check, K, p)

da <- simdataset(n, Pi = param$Pi, Mu = param$Mu, S = param$S)

save(da, file = "./data/simu.data.rda")
