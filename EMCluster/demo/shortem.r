library(EMCluster, quiet = TRUE)
set.seed(1234)

x2 <- da2$da
emobj <- list(pi = da2$pi, Mu = da2$Mu, LTSigma = da2$LTSigma)
ret.2 <- shortemcluster(x2, emobj, maxiter = 100, eps = 1e-2)
summary(ret.2)
