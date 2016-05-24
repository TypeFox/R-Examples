require(methods)
library(MixAll)

data(DebTrivedi)
dt <- DebTrivedi[1:500, c(1, 6,8, 15)]

model <- clusterPoisson( data=dt, nbCluster=2
                       , modelNames=clusterPoissonNames(prop = "equal")
                       , strategy = clusterFastStrategy())

model
