library(EMCluster, quiet = TRUE)

x2 <- da2$da
emobj <- list(pi = da2$pi, Mu = da2$Mu, LTSigma = da2$LTSigma)

eobj <- e.step(x2, emobj)
emobj <- m.step(x2, eobj)
summary(emobj)

eobj <- list(Gamma = class2Gamma(da2$class))
emobj <- m.step(x2, eobj)
summary(emobj)

