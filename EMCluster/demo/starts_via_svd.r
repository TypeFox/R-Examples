library(EMCluster, quiet = TRUE)
set.seed(1234)

nclass <- 10
x1 <- da1$da

ret.0 <- starts.via.svd(x1, nclass = nclass, method = "kmeans")
summary(ret.0)

ret.1 <- starts.via.svd(x1, nclass = nclass, method = "em")
summary(ret.1)

par(mfrow = c(1, 2))
plot2d(x1, color.pch = ret.0$class, main = ret.0$method)
plot2d(x1, color.pch = ret.1$class, main = ret.1$method)

