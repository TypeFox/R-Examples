library(parallel)
library(simPopulation)
library(party)
library(LiblineaR)
library(stringr)
library(microbenchmark)
library(e1071)

rm(list=ls())
seed <- 1234
data(eusilcS)   # load sample data

# 0.1s
system.time({
  eusilcP <- simStructure(eusilcS)
})

sapply(list.files("R", full.names=TRUE), source)


# non-parallel ~ 24secs
system.time({
  eusilcP2 <- simCategoricalOld(eusilcS, eusilcP)
})

# conditional probabilities
system.time({
  eusilcP2_0 <- simCategorical(eusilcS, eusilcP, method="distribution", parallel=TRUE)
})
# multinomial regression from package multinom
system.time({
  eusilcP2_1 <- simCategorical(eusilcS, eusilcP, method="multinom", parallel=TRUE)
})
# recursive partitioning from package party
#system.time({
#  eusilcP2_2 <- simCategorical(eusilcS, eusilcP, method="ctree", parallel=FALSE)
#})
# naivebayes from package e1071
system.time({
  eusilcP2_3 <- simCategorical(eusilcS, eusilcP, method="naivebayes", parallel=TRUE)
})
# liblinear from package LiblineaR
#system.time({
#  eusilcP2_4 <- simCategorical(eusilcS, eusilcP, method="liblinear", parallel=FALSE)
#})

# results
tab <- rbind(
  table(eusilcP2_0$pl030), # distr
  table(eusilcP2_1$pl030), # multinom
  #table(eusilcP2_2$pl030), # ctree
  table(eusilcP2_3$pl030) # nbayes
  #table(eusilcP2_4$pl030) # liblinear
)
barplot(tab, beside=TRUE, legend.text=c("distr","multinom","nbayes"), args.legend=list(x="topright"))

# benchmarking
nrruns <- 5
microbenchmark(
  simCategorical(eusilcS, eusilcP, method="distribution", parallel=TRUE),
  simCategorical(eusilcS, eusilcP, method="multinom", parallel=TRUE),
  #simCategorical(eusilcS, eusilcP, method="ctree", parallel=TRUE),
  simCategorical(eusilcS, eusilcP, method="naivebayes", parallel=TRUE),
  #simCategorical(eusilcS, eusilcP, method="liblinear", parallel=TRUE),
  times = 5
)


