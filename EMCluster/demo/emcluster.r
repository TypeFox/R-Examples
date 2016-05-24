library(EMCluster, quiet = TRUE)
set.seed(1234)

x2 <- da2$da
ret2.g <- emgroup(x2, nclass = 5)
ret2.gs <- summary(ret2.g)

emobj <- simple.init(x2, nclass = 5)
ret2.f <- emcluster(x2, emobj, assign.class = TRUE)
ret2.fs <- summary(ret2.f)

x3 <- da3$da
ret3.g <- emgroup(x3, nclass = 5)
ret3.gs <- summary(ret3.g)

emobj <- simple.init(x3, nclass = 5)
ret3.f <- emcluster(x3, emobj, assign.class = TRUE)
ret3.fs <- summary(ret3.f)

par(mfrow = c(2, 2))
plotem(ret2.gs, da2$da, main = "da2")
plotem(ret3.gs, da3$da, main = "da3")
plotem(ret2.fs, da2$da, main = "da2 simple")
plotem(ret3.fs, da3$da, main = "da3 simple")

