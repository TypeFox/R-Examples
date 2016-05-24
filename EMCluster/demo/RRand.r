library(EMCluster, quiet = TRUE)
set.seed(1234)

x2 <- da2$da
emobj <- emgroup(x2, nclass = 10)
ret.2 <- emcluster(x2, emobj, assign.class = TRUE)

RRand(da2$class, ret.2$class)

